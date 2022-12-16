from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

adenine = ["N9","C8","N7","C5","C4","C6","N3","C2","N1","N6","C1'"]
thymine = ["N1","C2","O2","N3","C4","O4","C5","C7","C6","C1'"]

delta_atoms = {"C5'": 0, "C4'": 1, "C3'": 2, "O3'": 3}

try:
  pdb = PDBFile(sys.argv[1])
except:
  raise

#forcefield = ForceField('amber14/DNA.OL15.xml', 'amber10_obc.xml')
forcefield = ForceField('amber14/DNA.OL15.xml')

#modeller = Modeller(pdb.topology, pdb.positions)
#modeller.addSolvent(forcefield, padding=1.0*nanometers, ionicStrength=0.1*molar, positiveIon='Na+')

# no constraints possible for atoms with 0 mass
#system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer, constraints=HBonds)

# explicit (vacuum)
system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer)
# implicit (no DNA params)
#system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer, implicitSolvent=OBC2, implicitSolventKappa=1.0/nanometer)

restraint = PeriodicTorsionForce()
system.addForce(restraint)

#integrator = LangevinMiddleIntegrator(50*kelvin, 1/picosecond, 0.001*picoseconds)
integrator = LangevinIntegrator(50*kelvin, 1/picosecond, 0.001*picoseconds)

try:
  platform = Platform.getPlatformByName(sys.argv[2])
  print("using", sys.argv[2], "platform")
except:
  platform = Platform.getPlatformByName('OpenCL')
  print("using OpenCL platform")

deltas = {}

print("Freezing base heavy atoms ...")
for atom in pdb.topology.atoms():
  if atom.name in delta_atoms:
    if atom.residue.chain.id+"_"+atom.residue.name+"_"+str(atom.residue.index) not in deltas:
      deltas[atom.residue.chain.id+"_"+atom.residue.name+"_"+str(atom.residue.index)] = [None, None, None, None]
    deltas[atom.residue.chain.id+"_"+atom.residue.name+"_"+str(atom.residue.index)][delta_atoms[atom.name]] = atom.index

  if atom.residue.name == "DA" and atom.name in adenine:
    #print("freezing atom " + atom.name)
    system.setParticleMass(atom.index,0.0)

  if atom.residue.name == "DT" and atom.name in thymine:
    system.setParticleMass(atom.index,0.0)

rest_count = 0

print("Setting delta restaints ...")
for delta in deltas:
  #print(delta, deltas[delta])
  if (deltas[delta][0] is not None and deltas[delta][1] is not None and deltas[delta][2] is not None and deltas[delta][3] is not None):
    restraint.addTorsion(deltas[delta][0], deltas[delta][1], deltas[delta][2], deltas[delta][3], 1, -0.7*radians, 1000*kilojoules_per_mole)
    rest_count += 1
  else:
    print("missing atoms in", delta)

print("Assigned", rest_count, "restraints")

simulation = Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

#simulation.reporters.append(PDBReporter('output.pdb', 10))
#simulation.reporters.append(StateDataReporter(stdout, 10, step=True, potentialEnergy=False, temperature=False))

ncycles=5

for cycle in range(ncycles):
  print('Running 20 step minimization batch ' + '%d' % (cycle+1) + " of " + str(ncycles) + " ...")

  simulation.minimizeEnergy(maxIterations=20, tolerance=1*kilojoule/mole)
  state = simulation.context.getState(getPositions=True, getEnergy=True, getForces=True)
  positions = state.getPositions()
  # (un)comment output of geometry snapshots from each batch
  #PDBFile.writeFile(simulation.topology, positions, open('output_' + '%02d' % cycle + '_min.pdb', 'w'))

PDBFile.writeFile(simulation.topology, positions, open('optimized_' + sys.argv[1], 'w'))

#  simulation.step(5)

#  state = simulation.context.getState(getPositions=True, getEnergy=True, getForces=True)
#  positions = state.getPositions()
#  PDBFile.writeFile(simulation.topology, positions, open('output_' + '%02d' % cycle + '_MD.pdb', 'w'))
