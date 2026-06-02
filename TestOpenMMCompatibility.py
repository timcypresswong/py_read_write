from openmm.app import PDBFile, ForceField, Simulation
from openmm import CustomExternalForce, Platform, VerletIntegrator, LangevinIntegrator
from openmm import unit as openmm_unit
 
import openmm

import argparse




def parser_args():
    parser = argparse.ArgumentParser(description='This is a simple program to test whether a pdb file can be run in openmm')
	
    parser.add_argument('pdb', help='input pdb name')
    parser.add_argument('--add_xml', help='additional xml file')

    args = parser.parse_args()
    return args

def run_openmm(pdb, add_xml):
    complex = PDBFile(pdb)
    force_field = ForceField('amber14-all.xml')
    if add_xml is not None:
        # force_field = ForceField( add_xml)
        # force_field = ForceField('amber14-all.xml', add_xml)
        force_field = ForceField('amber/protein.ff14SB.xml', add_xml)

    system = force_field.createSystem(
        complex.topology,
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
    simulation = Simulation(complex.topology,system, integrator=integrator, platform = platform)
    simulation.context.setPositions(complex.positions)
    simulation.minimizeEnergy(maxIterations = 10000)
    # better minimze
    minimized_positions = simulation.context.getState(getPositions = True).getPositions()

def run():
    args = parser_args()
    pdb = args.pdb
    add_xml = args.add_xml
    run_openmm(pdb, add_xml)

if __name__ == '__main__':
	run()