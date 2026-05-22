# from openmmforcefields.generators import GAFFTemplateGenerator
from openff.toolkit import Molecule
from openff.toolkit import ForceField, Molecule, Topology
from openff.interchange import Interchange
import argparse
import pickle

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
    generate_ff(mol,  base_name = base_name)

def generate_ff(mol, base_name):
    force_field = ForceField("openff-2.0.0.offxml")
    topology = mol.to_topology()
    interchange = Interchange.from_smirnoff(force_field = force_field, topology = topology)
    interchange.positions = mol.conformers[0]
    with open(base_name + "_interchange.pkl", 'wb') as f:
        pickle.dump(interchange, f)
    with open(base_name + ".offxml", "w") as f:
        f.write(force_field.to_string())


if __name__ == '__main__':
	run()