import argparse

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdEHTTools
from openmm.app import PDBFile, Simulation, Modeller
from openff.toolkit import ForceField, Topology, Molecule
from openmm import CustomExternalForce, Platform, VerletIntegrator, LangevinIntegrator
from openmm import unit as openmm_unit
import openmm
from openmm import Vec3
from openff.toolkit import Molecule as OFFMolecule
import pickle
import GensmallFF_XML

def parser_args():
    parser = argparse.ArgumentParser(description='Join the protein and ligand to a complex pdb file')
	
    parser.add_argument('protein', help='input protein pdb name')
    parser.add_argument('ligand', help='input ligand mol/sdf name')
    parser.add_argument('complex', help='output complex pdb name')
    args = parser.parse_args()
    return args

def join_protein_ligand_fast(protein_path, ligand_path, complex_path): # this generate forcefield for the entire complex using interchange method, for a large complex system, it would be very slow
    protein = PDBFile(protein_path)
    ligand = Chem.MolFromMolFile(ligand_path, removeHs = False, sanitize = True)
    off_ligand = OFFMolecule.from_rdkit(ligand)

    ligand_top = off_ligand.to_topology().to_openmm()
    ligand_RDpos = off_ligand.conformers[0]

    ligand_pos = [  openmm_unit.Quantity(Vec3( x,y,z ), openmm_unit.angstrom) for x,y,z in ligand_RDpos.magnitude  ]


    complex_model = Modeller(protein.topology, protein.positions)
    complex_model.add(ligand_top, ligand_pos)

    with open(complex_path, 'w') as f:
        PDBFile.writeFile(complex_model.topology, complex_model.positions, f)

    base_name = complex_path.split('_')[0]
    GensmallFF_XML.parmed_generate_ff(off_ligand, base_name + "_drug")

    force_field = openmm.app.ForceField('amber/protein.ff14SB.xml', base_name + "_drug_openff.xml")
    system = force_field.createSystem(
        complex_model.topology,
        constraints = None,
    )
    integrator = openmm.LangevinIntegrator(
    300 * openmm_unit.kelvin,
    1 / openmm_unit.picosecond,
    0.002 * openmm_unit.picoseconds,
    )
    try:
        platform = Platform.getPlatformByName("CUDA")
    except:
        platform = Platform.getPlatformByName("CPU")
    simulation = Simulation(complex_model.topology,system, integrator=integrator, platform = platform)
    simulation.context.setPositions(complex_model.positions)
    simulation.minimizeEnergy(maxIterations = 10000)

    minimized_positions = simulation.context.getState(getPositions = True).getPositions()

    with open(base_name + "off.pdb", 'w') as f:
        PDBFile.writeFile(simulation.topology, minimized_positions, f)


def join_protein_ligand_slow(protein_path, ligand_path, complex_path):
    protein = PDBFile(protein_path)
    ligand = Chem.MolFromMolFile(ligand_path, removeHs = False, sanitize = True)
    off_ligand = OFFMolecule.from_rdkit(ligand)

    ligand_top = off_ligand.to_topology().to_openmm()
    ligand_RDpos = off_ligand.conformers[0]
    # print(type(protein.positions))
    # print(type(protein.positions[0]))
    # print(type(ligand_RDpos))
    # print(type(ligand_RDpos[0]))
    # print(type(ligand_pos))

    # print(ligand_RDpos.magnitude)
    # print(ligand_RDpos.units)

    # print(protein.positions[0])
    # print(ligand_pos[0])



    ligand_pos = [  openmm_unit.Quantity(Vec3( x,y,z ), openmm_unit.angstrom) for x,y,z in ligand_RDpos.magnitude  ]
    complex_model = Modeller(protein.topology, protein.positions)
    complex_model.add(ligand_top, ligand_pos)




    # check accidental bond
    # protein_atoms = list( protein.topology.atoms() )
    # num_atom_protein = len(protein_atoms)
    # for bond in complex_model.topology.bonds():
    #     idx1 = bond[0].index
    #     idx2 = bond[1].index
    #     if (idx1 < num_atom_protein and idx2 >= num_atom_protein) or (idx2 < num_atom_protein and idx1 >= num_atom_protein):
    #         print("found accidental bond!", idx1, idx2)

    # force_field = ForceField("openff-2.0.0.offxml", 'ff14sb_off_impropers_0.0.2.offxml', 'ff14sb_0.0.2.offxml', 'GBSA_OBC2-1.0.offxml')

    force_field = ForceField("openff-2.2.0.offxml", 'ff14sb_off_impropers_0.0.4.offxml')

    OFF_Protein_top = Topology.from_pdb(protein_path)
    OFF_ligand = Molecule.from_file(ligand_path)
    new_top_mols = []
    for molecule in OFF_Protein_top.molecules:
        new_top_mols.append(molecule)
    new_top_mols.append(OFF_ligand)
    OFF_Complex_top = Topology.from_molecules(new_top_mols)

    interchange = force_field.create_interchange( topology = OFF_Complex_top)

    base_name = complex_path.rsplit('.', maxsplit = 1)[0]
    # with open(base_name + "interchange.pkl", 'wb') as f:
    #     pickle.dump(interchange, f)

    with open(complex_path, 'w') as f:
        PDBFile.writeFile(complex_model.topology, complex_model.positions, f)

    integrator = openmm.LangevinIntegrator(
    300 * openmm_unit.kelvin,
    1 / openmm_unit.picosecond,
    0.002 * openmm_unit.picoseconds,
    )
    try:
        platform = Platform.getPlatformByName("CUDA")
    except:
        platform = Platform.getPlatformByName("CPU")
    simulation = interchange.to_openmm_simulation(integrator=integrator, platform = platform)
    simulation.context.setPositions(complex_model.positions)
    simulation.minimizeEnergy(maxIterations = 10000)
    # better minimze
    minimized_positions = simulation.context.getState(getPositions = True).getPositions()

    with open(base_name + "off.pdb", 'w') as f:
        PDBFile.writeFile(simulation.topology, minimized_positions, f)


def run():
    args = parser_args()
    protein_path = args.protein
    ligand_path = args.ligand
    complex_path = args.complex

    join_protein_ligand_fast(protein_path, ligand_path, complex_path)

if __name__ == '__main__':
	run()

    