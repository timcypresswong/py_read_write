# from openmmforcefields.generators import GAFFTemplateGenerator
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
import openff.toolkit
from openff.toolkit import  Molecule, Topology #ForceField,
from openff.interchange import Interchange
from openmm import app, XmlSerializer
from openff.units.openmm import ensure_quantity
from openff.units import unit
import openmm.app
import parmed
import argparse
import pickle
import xml.etree.ElementTree as ET
from importlib import resources
import pathlib

def parser_args():
    parser = argparse.ArgumentParser(description='This utilzize openff to generate .xml forcefields for small molecules. Currently only support openff and gaff forcefield')
    parser.add_argument('filename', help='input mol or sdf name')
    args = parser.parse_args()
    return args

def run():
    args = parser_args()
    filename = args.filename
    base_name = filename.rsplit('.', maxsplit = 1)[0]
    mol = Molecule.from_file(filename)
    # generate_ff(mol,  base_name = base_name)
    parmed_generate_ff(mol,  base_name = base_name)

def generate_ff(mol, base_name):
    force_field = openff.toolkit.ForceField("openff-2.0.0.offxml")
    topology = mol.to_topology()
    interchange = Interchange.from_smirnoff(force_field = force_field, topology = topology)
    interchange.positions = mol.conformers[0]
    # with open(base_name + "_interchange.pkl", 'wb') as f:
    #     pickle.dump(interchange, f)
    with open(base_name + ".offxml", "w") as f:
        f.write(force_field.to_string())

def parmed_generate_ff(mol, base_name):
    smirnoff = SMIRNOFFTemplateGenerator(
        molecules = mol,
        forcefield = "openff-2.3.0.offxml"
    )

    # tempalte_file_name = 'amber14-all.xml'
    tempalte_file_name = "amber14/protein.ff14SB.xml"

    # force_field = openmm.app.ForceField()
    force_field = openmm.app.ForceField(tempalte_file_name)
    # force_field = openmm.app.ForceField('amber/protein.ff14SB.xml')

    force_field.registerTemplateGenerator( smirnoff.generator )

    topology = mol.to_topology()
    topology = topology.to_openmm()
    system = force_field.createSystem(
        topology,
        nonbondedCutoff = ensure_quantity(0.9 * unit.nanometer, "openmm"),
        switchDistance = ensure_quantity(0.8 * unit.nanometer, "openmm")
        )

    # nonbonded_force = None
    # for i in range(system.getNumForces()):
    #     if isinstance(system.getForce(i), openmm.NonbondedForce):
    #         nonbonded_force = system.getForce(i)
    #         break
    # if nonbonded_force is not None:
    #     # 2. 设置你需要的缩放因子
    #     # 设定 Coulomb 14 scaling factor = 0.833333
    #     nonbonded_force.setCoulomb14Scale(0.833333)
    #     # 设定 LJ 14 scaling factor = 0.5
    #     nonbonded_force.setLJ14Scale(0.5)

    st = parmed.openmm.load_topology(topology, system = system)
    w = parmed.amber.parameters.ParameterSet.from_structure(st)
    ww = parmed.openmm.parameters.OpenMMParameterSet.from_parameterset(w)
    ww.residues.update(parmed.modeller.ResidueTemplateContainer.from_structure(st).to_library())
    ww.write(base_name + "_openff.xml")
    modified_LJCoulomb_scaling(tempalte_file_name, base_name + "_openff.xml")

def modified_LJCoulomb_scaling(tempalte_file_name, to_modified_filename):


    app_pkg_path = pathlib.Path(openmm.app.__file__).parent
    data_dir = app_pkg_path / 'data'  # 内置力场文件目录

    # 定位目标文件
    ff_path = data_dir / tempalte_file_name
    # print(ff_path)

    tree = ET.parse(ff_path)
    root = tree.getroot()
    nonbonded_elem = root.find('.//NonbondedForce')
    coul14 = None
    lj14 = None
    if nonbonded_elem is not None:
        coul14 = nonbonded_elem.get('coulomb14scale')
        lj14 = nonbonded_elem.get('lj14scale')
        # print(f"coulomb14scale = {coul14}, lj14scale = {lj14}")
    else:
        print("在 XML 中未找到 <NonbondedForce> 元素。")


    tree = ET.parse(to_modified_filename)
    root = tree.getroot()

    # 查找 NonbondedForce 元素并修改其属性
    for nonbonded in root.findall('.//NonbondedForce'):
        if coul14 is not None:
            nonbonded.set('coulomb14scale', coul14)
        if lj14 is not None:
            nonbonded.set('lj14scale', lj14)

    # 将修改后的内容写回文件
    tree.write(to_modified_filename, encoding='UTF-8', xml_declaration=True)



if __name__ == '__main__':
	run()