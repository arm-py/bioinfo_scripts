from openmm.app import *
from openmm import *
from openmm.unit import *
from pdbfixer import PDBFixer # Importamos el arreglador
from sys import stdout

# 1. Cargamos y arreglamos el archivo
print("Arreglando el archivo PDB (añadiendo hidrógenos)...")
fixer = PDBFixer('receptor-ligand_model1.pdb')
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0) # Añade hidrógenos para pH neutro

# 2. Configuración del sistema
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
# Usamos la topología y posiciones del 'fixer'
system = forcefield.createSystem(fixer.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(fixer.topology, system, integrator)
simulation.context.setPositions(fixer.positions)

# 3. Minimización y ejecución
print("Minimizando energía...")
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

# --- LÍNEA NUEVA AÑADIDA AQUÍ ---
simulation.reporters.append(DCDReporter('trayectoria.dcd', 1000))
# --------------------------------

print("Iniciando dinámica...")
simulation.step(10000)
print("Simulación terminada. Revisa el archivo output.pdb y trayectoria.dcd")