import opengate as gate
from opengate.tests import utility
from scipy.spatial.transform import Rotation as R
import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import pyvista as pv
from functions import DD_prob
import spekpy as sp

sim = gate.Simulation()
sim.volume_manager.add_material_database("/home/ldauriac/Software/opengate/opengate/opengate/tests/data/GateMaterials.db")

# units
m = gate.g4_units.m
cm = gate.g4_units.cm
mm = gate.g4_units.mm
um = gate.g4_units.um
keV = gate.g4_units.keV
MeV = gate.g4_units.MeV
Bq = gate.g4_units.Bq

world = sim.world
world.size = [2*m, 2*m, 2*m]
world.material = "G4_Galactic"

#volume correspondant au microarray contenant les samples
membrane = sim.add_volume("BoxVolume" , "membrane")
membrane.size = [26*mm, 76*mm , 0.00014*mm] 
membrane.translation = [0, 0, 0]
membrane.material = "G4_POLYETHYLENE" 
rotation_matrix_obj = R.from_euler('y', 15, degrees=True).as_matrix()
membrane.rotation = rotation_matrix_obj

cell_vol = sim.add_volume("EllipsoidVolume" , "cell_vol") #volume correspondant à la cellule collée à un écran de Polyéthylene
cell_vol.a = 26*mm
cell_vol.b = 76*mm
cell_vol.c = 0.1*mm
cell_vol.material = "G4_WATER"
cell_vol.mother = membrane
cell_vol.translation = [0,0, 0.05007*mm]

sample_témoin = sim.add_volume("TubsVolume" , "sample")
sample_témoin.rmin = 0
sample_témoin.rmax = 0.005*mm
sample_témoin.dz = 0.02*mm
sample_témoin.mother = cell_vol
sample_témoin.translation = [0,1*mm,0]
sample_témoin.material = "G4_WATER"
rotation_matrix_sample = R.from_euler('z', 90, degrees=True).as_matrix()
sample_témoin.rotation = rotation_matrix_sample

sample = sim.add_volume("TubsVolume" , "sample")
sample.rmin = 0
sample.rmax = 0.005*mm
sample.dz = 0.02*mm
sample.mother = cell_vol
sample.translation = [0,0,0]
sample.material = "G4_Water1ppmPt"
rotation_matrix_sample = R.from_euler('z', 90, degrees=True).as_matrix()
sample.rotation = rotation_matrix_sample

sample2 = sim.add_volume("TubsVolume" , "sample")
sample2.rmin = 0
sample2.rmax = 0.005*mm
sample2.dz = 0.02*mm
sample2.mother = cell_vol
sample2.translation = [0,-1*mm,0]
sample2.material = "G4_3ppmPt"
rotation_matrix_sample = R.from_euler('z', 90, degrees=True).as_matrix()
sample2.rotation = rotation_matrix_sample

#sphère de détection en Ge pour faciliter la mesure
spheric_detector = sim.add_volume("SphereVolume", "spheric_detector")
spheric_detector.mother = "vacoat"
spheric_detector.rmin = 38.99*mm
spheric_detector.rmax = 39*mm
spheric_detector.material = "G4_Ge"

#détecteur HPGe
HPGe = sim.add_volume("TubsVolume" , "HPGe")
HPGe.rmin = 0
HPGe.rmax = 5*mm
HPGe.dz = 3*mm
HPGe.material = "G4_Ge"
HPGe.translation = [-7.95*cm, 0, 0] #6.5cm+13mm+1.5mm dans la manip 
rotation_matrix_detector = R.from_euler('y', 90, degrees=True).as_matrix()
HPGe.rotation = rotation_matrix_detector #rotation de 90° sur l'axe X pour que la surface sensible du détecteur soit face au microarray.

#actor
HPGe_detection = sim.add_actor("PhaseSpaceActor", "HPGE_spheric_detection")
HPGe_detection.attached_to = ["spheric_detector"]
HPGe_detection.attributes = ['KineticEnergy']
HPGe_detection.output_filename = 'HPGE__spheric_detection.root'

#filter
parfilter = sim.add_filter("ParticleFilter", "parfilter")
parfilter.particle = "gamma"
HPGe_detection.filters.append(parfilter)


#physique
sim.physics_manager.physics_list_name = "G4EmStandardPhysics_option4"
sim.physics_manager.em_parameters.fluo = True
sim.physics_manager.em_parameters.auger = True  
sim.physics_manager.em_parameters.auger_cascade = True
sim.physics_manager.global_production_cuts.gamma = 1 * um
sim.physics_manager.global_production_cuts.electron = 100 * um

#source
"""spek = sp.Spek()
spek.set_kvp = 100
spek.set_bw = 0.35
spek.set_targ ="W" 
spectrum = spek.get_spectrum()"""

x_positions = np.linspace(-3, 3, 4) 
y_positions = np.linspace(-3, 3, 4)
translations = [[x, y, -0.5] for x in x_positions for y in y_positions]

for translation in translations : #exemples de translation 
  beam_name = f"beam_{translation[0]}_{translation[1]}_{translation[2]}" 
  incident_beam = sim.add_source('GenericSource' , beam_name)
  """incident_beam.energy.type = "spectrum_discrete"
  incident_beam.energy.spectrum_energies = spek.get_spectrum()[0] #énergies possibles du faisceau poly
  incident_beam.energy.spectrum_weights = spek.get_spectrum()[1]  #poids des énergies du faisceau poly"""
  #incident_beam.position.type = "box"
  #incident_beam.position.half_size = [10*mm,10*mm, 0]
  incident_beam.particle = "gamma"
  incident_beam.position.type = "point"
  incident_beam.direction.type = "momentum"
  incident_beam.direction.momentum = [0, 0, 1]
  incident_beam.attached_to = "sample" 
  incident_beam.position.translation = translation
  incident_beam.energy.mono = 100 * keV
  incident_beam.n = 100


stats = sim.add_actor("SimulationStatisticsActor","stats")
stats.track_types_flag = True

sim.visu = True
sim.visu_type = "vrml_file_only"
sim.visu_filename = "simulation.wrl"
sim.visu_verbose = False

# === Exécuter la simulation pour générer le fichier VRML ===
print("🔄 Exécution de la simulation pour générer la visualisation...")
sim.run()
print("✅ Simulation terminée.")

# === Charger le fichier VRML avec PyVista ===
print("📂 Chargement du fichier VRML pour affichage...")
pl = pv.Plotter()
pl.import_vrml("simulation.wrl")  # Charger la scène OpenGATE

# === Ajouter des axes colorés pour l'orientation ===
axes = pv.Axes()
axes.axes_actor.total_length = [1, 1, 1]  # Taille des axes
axes.axes_actor.x_axis_shaft_properties.color = (1, 0, 0)  # Rouge pour X
axes.axes_actor.y_axis_shaft_properties.color = (0, 1, 0)  # Vert pour Y
axes.axes_actor.z_axis_shaft_properties.color = (0, 0, 1)  # Bleu pour Z
pl.add_actor(axes.axes_actor)

# === Paramètres d'affichage ===
pl.background_color = "black"  # Fond noir pour meilleur contraste
pl.show_grid()  # Afficher une grille
pl.show()  # Afficher la scène 3D

print("lancement de la simulation")
print(sim.volume_manager.volumes.keys())  # Liste tous les volumes définis
detection_prob = DD_prob(incident_beam.n)
prob_tot = detection_prob[0]
detected_particles = detection_prob[1]
print(f"la probabilité de détection d'une particule est {prob_tot}")
print(f"Pour {incident_beam.n} on doit détecter {detected_particles}")

#sim.run()
print("fin de la simulation")
print(stats)
print(f"nombre total de particules primaires émises:", stats.counts.events)