#initial 200 ns simulations (2 fs timestep, 1 nm non-bonded cutoff)

from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import time

start_time = time.time()

prmtop = AmberPrmtopFile('input.prmtop')
inpcrd = AmberInpcrdFile('input.inpcrd')

system = prmtop.createSystem(nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=1*nanometer,
        constraints=HBonds, implicitSolvent=GBn2)
integrator = LangevinMiddleIntegrator(298*kelvin, 1/picosecond, 0.002*picoseconds)
platform = Platform.getPlatformByName('CUDA')

simulation = Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(inpcrd.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(100000000)

print ("simulation time:", time.time() - start_time, "s")
