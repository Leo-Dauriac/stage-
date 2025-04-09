import opengate as gate
#from opengate.tests import utility
#from scipy.spatial.transform import Rotation as R
import numpy as np 
import matplotlib.pyplot as plt 
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import pyvista as pv
from functions import choix_element, is_overlapping, random_point_in_cylinder, DD_prob, test_raie_prob
import spekpy as sp


sim = gate.Simulation()

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

sphere_Ag = sim.add_volume("SphereVolume", "sphere")
sphere_Ag.rmin = 0
sphere_Ag.rmax = 1*mm
sphere_Ag.material = "G4_WATER"
sphere_Ag.translation = [0, 0, 0]

#sphère de détection en Ge pour faciliter la mesure
spheric_detector = sim.add_volume("SphereVolume", "spheric_detector")
spheric_detector.mother = "world"
spheric_detector.rmin = 1.9*mm
spheric_detector.rmax = 2*mm
spheric_detector.material = "G4_WATER"

spheric_hits = sim.add_actor('DigitizerHitsCollectionActor', 'test_raie_hits')
spheric_hits.attached_to = "spheric_detector"
spheric_hits.output_filename = 'test_raie_hits.root'
spheric_hits.attributes = ['TotalEnergyDeposit']

ps_actor = sim.add_actor("PhaseSpaceActor", "PS_out")
ps_actor.attached_to = "spheric_detector"  # Enregistre tout ce qui sort de la sphère
ps_actor.output_filename = "fluorescence_out.root"
ps_actor.attributes = ["KineticEnergy"]

sim.physics_manager.physics_list_name = "G4EmStandardPhysics_option4"
sim.physics_manager.em_parameters.fluo = True
sim.physics_manager.em_parameters.auger = True  
sim.physics_manager.em_parameters.auger_cascade = True
sim.physics_manager.global_production_cuts.gamma = 0.001 * um
sim.physics_manager.global_production_cuts.electron = 100 * um

"""spek = sp.Spek()
spek.set_kvp = 100
spek.set_bw = 0.35
spek.targ ="W"""
"""energies = spek.get_spectrum()[0]
intensities = spek.get_spectrum()[1]
print(spectrum)
plt.figure(figsize=(8, 5))
plt.plot(energies, intensities, label="Spectre X 100 kV - Tungstène", color="blue")
# Personnalisation du graphique
plt.xlabel("Énergie des photons (keV)")
plt.ylabel("Intensité relative")
plt.title("Spectre X généré avec SpekPy")
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()
# Affichage
plt.show()"""
x_positions = np.linspace(-0.5, 0.5, 4) 
y_positions = np.linspace(-0.5, 0.5, 4)
translations = [[x, y, -0.12] for x in x_positions for y in y_positions]

for translation in translations : #exemples de translation 
  beam_name = f"beam_{translation[0]}_{translation[1]}_{translation[2]}" 
  incident_beam = sim.add_source('GenericSource' , beam_name)
#incident_beam = sim.add_source('GenericSource' , 'MyFirstSource')
  incident_beam.position.type = "point"
  incident_beam.particle = "gamma"
#incident_beam.position.half_size = [0.5*mm,0.2*mm, 0]
  incident_beam.direction.type = "momentum"
  incident_beam.direction.momentum = [0, 0, 1]
  incident_beam.attached_to = "spheric_detector"
  #incident_beam.position.translation = [0*mm, 0*mm, -0.12*mm]
  incident_beam.position.translation = translation
  incident_beam.n = 100000
  incident_beam.energy.mono = 69 * keV



stats = sim.add_actor("SimulationStatisticsActor","stats")
stats.track_types_flag = True

sim.visu = True
sim.visu_type = "vrml_file_only"
sim.visu_filename = "test.wrl"



# === Exécuter la simulation pour générer le fichier VRML ===
"""print("🔄 Exécution de la simulation pour générer la visualisation...")
sim.run()
print("✅ Simulation terminée.")

# === Charger le fichier VRML avec PyVista ===
print("📂 Chargement du fichier VRML pour affichage...")
pl = pv.Plotter()
pl.import_vrml("test.wrl")  # Charger la scène OpenGATE

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
pl.show()  # Afficher la scène 3D"""

detection_prob = test_raie_prob(incident_beam.n)
prob_tot = detection_prob[0]
detected_particles = detection_prob[1]
print(f"la probabilité de détection d'une particule est {prob_tot}")
print(f"Pour {incident_beam.n} on doit détecter {detected_particles}") #ne prendra pas en compte le nombre de photons fluo absorbé dans le volume d'Ag avanr d'être détecté
print("lancement de la simulation")
sim.run()
print("fin de la simulation")
print(stats)
print(f"nombre total de particules primaires émises:", stats.counts.events)

