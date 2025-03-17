import opengate as gate
from opengate.tests import utility
from scipy.spatial.transform import Rotation as R
import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import math
import xraylib 
from functions import choix_element, is_overlapping, random_point_in_cylinder, DD_prob



sim = gate.Simulation()
#paths = utility.get_default_test_paths(__file__, "gate_test009_voxels", "test020")
#sim.volume_manager.add_material_database(paths.data / "GateMaterials.db")

# units
m = gate.g4_units.m
cm = gate.g4_units.cm
mm = gate.g4_units.mm
keV = gate.g4_units.keV
MeV = gate.g4_units.MeV
Bq = gate.g4_units.Bq

#volume correspondant au microarray contenant les samples
microarray = sim.add_volume("BoxVolume" , "microarray")
microarray.size = [8*mm, 8*mm , 0.1*mm]
microarray.translation = [0, 0, 0]
microarray.material = "G4_PLEXIGLASS" 
rotation_matrix = R.from_euler('y', 15, degrees=True).as_matrix()
microarray.rotation = rotation_matrix #rotation de 15° autour de l'axe Y



#volume des échantillons cellulaires contenant les éléments fluorescents
#cela doit être un cylindre de 1mm de diamètre et 10µm de hauteur car la goutte est écrasée
i = 0 
j = 0
for i in range(8) : 
  for j in range(8) : 
    sample_name = f"sample_{i}_{j}" #pour donner un nom différent à chaque sample/goutelette
    sample = sim.add_volume("TubsVolume" , sample_name)
    sample.mother = "microarray"
    sample.material = "G4_WATER"
    sample.rmin = 0
    sample.rmax = 0.5*mm
    sample.dz = 5*1e-3*mm
    sample.translation = [
    (-3.5 + i)*1*mm ,
    (-3.5 + j)*1*mm ,
    0
     ]
    j+=1
  i+=1  

#source X-ray monoénergétique émetteur perpendiculaire au sample
#(translaté sur le premier pixel en ht à gche du sample, reste à le déplacer pour obtenir un balayage)
incident_beam = sim.add_source('GenericSource' , 'MyFirstSource')
incident_beam.position.type = "box"
incident_beam.position.half_size = [250*1e-6*mm,250*1e-6*mm, 0]
incident_beam.direction.type = "momentum"
incident_beam.direction.momentum = [0, 0, 1]
incident_beam.attached_to = "microarray"
incident_beam.position.translation = [0*mm, 0*mm, 0]
#incident_beam.activity = 100*Bq
incident_beam.n = 100000 
incident_beam.energy.mono = 69 * keV
#ou alors un faisceau type iso qui couvre 500X500nm mais dans ce cas plus de balayage

i  = 0
max_attempts = 1000  # Limite d'essais pour éviter un blocage
attempts = 0
existing_positions = []
for i in range(1000):
    z_tag_name = f"z_tag_{i}"
    z_tag = sim.add_volume("SphereVolume", z_tag_name)
    z_tag.rmin = 0
    z_tag.rmax = 0.005*mm
    parent_i = np.random.randint(0, 8)
    parent_j = np.random.randint(0, 8)
    z_tag.mother = f"sample_{parent_i}_{parent_j}" #permet d'attacher le z_tag à l'un des 64 samples aléatoirement
    
    while True:                               #génère de nouvelles coordonnées tant que is_overlapping est vraie
        new_pos = random_point_in_cylinder(0.5, 0.01, 0.005)

        # Vérifie si la position est valide
        if not is_overlapping(new_pos, existing_positions, 0.005):        
            #print(f"✅ Position acceptée après {attempts} essais : {new_pos}")
            break  # Position valide, on sort de la boucle
        attempts += 1

        if attempts == max_attempts:
         print("⚠️ Avertissement : Impossible de placer un `z_tag`, trop d'overlaps !")
         print(f"Tentative échouée après {attempts} essais.")
         print(f"Nombre de `z_tag` déjà placés : {len(existing_positions)}")

    # Ajouter la position validée à la liste
    existing_positions.append(new_pos)

    z_tag.translation = new_pos
    mat = choix_element()
    z_tag.material = mat
    i+=1



"""for i in range(10000) :
    z_tag_name = f"z_tag_{i}"
    z_tag = sim.add_volume("SphereVolume", z_tag_name)
    z_tag.rmin = 0
    z_tag.rmax = 5*10e-3*mm
    parent_i = np.random.randint(0, 8)
    parent_j = np.random.randint(0, 8)
    z_tag.mother = f"sample_{parent_i}_{parent_j}" #permet d'attacher le z_tag à l'un des 64 samples aléatoirement
    z_tag.translation = [
       ( np.random.uniform(0,10e-2)*mm ), 
       (np.random.uniform(0,10e-2)*mm),
       (np.random.uniform(0,10e-3)*mm)
    ]
    mat = choix_element()
    z_tag.material = mat
    i+=1"""

#print(sim.volume_manager.dump_volume_tree())    

#création des volumes détecteurs

"""#SDD1 au dessus du sample
SDD1 = sim.add_volume("TubsVolume" , "SDD1")  #tourner le cylindre pour que sa face circulaire soit vers le microarray ?
SDD1.mother = "world"
SDD1.rmin = 0
SDD1.rmax = 3.1*mm
SDD1.dz = 0.25 * mm
SDD1.material = "G4_Si"
SDD1.translation = [-34.25*mm, 0, 0] #l'isocentre du détecteur est à gauche du microarray sur l'axe x à une distance de 30mm 
#+0.25mm pour que la surface et non l'isocentre du détecteur soit à 30mm
#+4mm pour la moitié de la largeur du microarray 

#SDD2 au dessous du sample
SDD2 = sim.add_volume("TubsVolume" , "SDD2")
SDD2.mother = "world"
SDD2.rmin = 0
SDD2.rmax = 3.1*mm
SDD2.dz = 0.25 * mm
SDD2.material = "G4_Si"
SDD2.translation = [34.25*mm, 0, 0] #pareil mais détecteur à droite """

#sphère de détection en Ge pour faciliter la mesure
spheric_detector = sim.add_volume("SphereVolume", "spheric_detector")
spheric_detector.mother = "world"
spheric_detector.rmin = 34*mm
spheric_detector.rmax = 35*mm
spheric_detector.material = "G4_Ge"



"""print(f"Microarray position: {microarray.translation}")
print(f"SDD1 position: {SDD1.translation}, dz={SDD1.dz}")
print(f"SDD2 position: {SDD2.translation}, dz={SDD2.dz}")"""

"""#détecteur GeCMOS
GeCMOS = sim.add_volume("TubsVolume" , "GeCMOS")
GeCMOS.rmin = 0
GeCMOS.rmax = 5*mm
GeCMOS.dz = 0.5 * mm
GeCMOS.material = "G4_Ge"
GeCMOS.translation = [-34.5, 0, 0]"""

#Actors :

#Actor pour calculer la dose reçu par le microarray
dose_microarray = sim.add_actor("DoseActor" , "dose_microarray")
dose_microarray.attached_to = "microarray"
dose_microarray.size = [8*mm, 8*mm , 1*mm]
dose_microarray.spacing = [0.01*mm, 0.01*mm, 0.01*mm] #taille d'un pixel correspondant au microarray découpé en 8x8 pixels
#dose_microarray.dose = True
dose_microarray.output_filename = "microarray-dose.mhd"

hit_collec = sim.add_actor('DigitizerHitsCollectionActor', 'Hit_collec')
hit_collec.attached_to = "microarray"
hit_collec.output_filename = 'test_hits.root'
hit_collec.attributes = ['TotalEnergyDeposit', 'KineticEnergy', 'PostPosition','GlobalTime', 'RunID', 'ThreadID', 'TrackID']

spheric_hits = sim.add_actor('DigitizerHitsCollectionActor', 'spheric_detector_hits')
spheric_hits.attached_to = ["spheric_detector"]
spheric_hits.output_filename = 'spheric_detector_hits.root'
spheric_hits.attributes = ['TotalEnergyDeposit']

"""#energie déposée dans les deux détecteurs SDD
edep_detector = sim.add_actor("DoseActor", "edep-detector")
edep_detector.attached_to = ["SDD1", "SDD2"]
edep_detector.size = [6*mm, 6*mm, 1*mm]
edep_detector.spacing = [1*mm, 1*mm, 1*mm]
#edep_detector.edep_active = True
edep_detector.output_filename = "detector-edep.root"""

#energie déposée dans le deuxième détecteur SDD
#edep_detector2 = sim.add_actor("DoseActor", "edep-detector2")
#edep_detector2.attached_to = "SDD2"
#edep_detector2.rmin = 0
#edep_detector2.rmax = 3,1*mmexamples/extended/hadronic/Hadr09
#edep_detector2.height = 0,5*mm
#edep_detector2.edep_active = True
#edep_detector2.output_filename = "detector2-edep.mhd"

#LETACTOR
"""LET = sim.add_actor("LETActor" , "LET_sample")
LET.attached_to = ["microarray" , "SDD1" , "SDD2"]
LET.output_filename = "LET_output.mhd"""



#infos particules détecteur
"""phsp = sim.add_actor("PhaseSpaceActor", "phsp_detector")
phsp.attached_to = ["SDD1" , "SDD2"]
phsp.attributes = [
    "KineticEnergy", #in fine isoler le pic PE des éléments
    "Weight", #analyser la contribution du diffusé ?
    "PostPosition", #preposition extérieure/postposition intérieure = particule entrante
    "PrePosition", #preposition intérieure/postposition extérieure = particule sortante
    "ParticleName", 
    "PreDirection",
    "PostDirection",
    "EventPosition", #contribution de la source et du fond
]
phsp.output_filename = "phsp_detector.root"""

stats = sim.add_actor("SimulationStatisticsActor","stats")
stats.track_types_flag = True

detection_prob = DD_prob(incident_beam.n)
prob_tot = detection_prob[0]
detected_particles = detection_prob[1]
print(f"la probabilité de détection d'une particule est {prob_tot}")
print(f"Pour {incident_beam.n} on doit détecter {detected_particles}")

print("lancement de la simulation")
#print(sim.volume_manager.volumes.keys())  # Liste tous les volumes définis

sim.run()
print("fin de la simulation")
print(stats)
print(f"nombre total de particules primaires émises:", stats.counts.events)





