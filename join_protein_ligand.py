import argparse

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdEHTTools
from openmm.app import PDBFile, ForceField, Simulation, Modeller
from openmm import CustomExternalForce, Platform, VerletIntegrator, LangevinIntegrator
from openmm import unit as openmm_unit
import openmm
from openmm import Vec3
from openff.toolkit import Molecule as OFFMolecule


def parser_args():
    parser = argparse.ArgumentParser(description='Join the protein and ligand to a complex pdb file')
	
    parser.add_argument('protein', help='input protein pdb name')
    parser.add_argument('ligand', help='input ligand mol/sdf name')
    parser.add_argument('complex', help='output complex pdb name')
    args = parser.parse_args()
    return args

def join_protein_ligand(protein_path, ligand_path, complex_path):
    protein = PDBFile(protein_path)
    ligand = Chem.MolFromMolFile(ligand_path, removeHs = False, sanitize = True)
    off_ligand = OFFMolecule.from_rdkit(ligand)

    ligand_top = off_ligand.to_topology().to_openmm()
    ligand_RDpos = off_ligand.conformers[0]

    ligand_pos = [  openmm_unit.Quantity(Vec3( x,y,z ), openmm_unit.angstrom) for x,y,z in ligand_RDpos.magnitude  ]

    # print(type(protein.positions))
    # print(type(protein.positions[0]))
    # print(type(ligand_RDpos))
    # print(type(ligand_RDpos[0]))
    # print(type(ligand_pos))

    # print(ligand_RDpos.magnitude)
    # print(ligand_RDpos.units)

    # print(protein.positions[0])
    # print(ligand_pos[0])


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

    with open(complex_path, 'w') as f:
        PDBFile.writeFile(complex_model.topology, complex_model.positions, f)



def run():
    args = parser_args()
    protein_path = args.protein
    ligand_path = args.ligand
    complex_path = args.complex

    join_protein_ligand(protein_path, ligand_path, complex_path)

if __name__ == '__main__':
	run()

    