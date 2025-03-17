import opengate as gate
from opengate.tests import utility
from scipy.spatial.transform import Rotation
import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import math
import uproot

sim = gate.Simulation()


# Dimensions du tableau (correspondant à DimSize dans le fichier .mhd)
dim_x, dim_y, dim_z = 8, 8, 1

# Chargement des données en tant que tableau NumPy
data = np.fromfile("microarray-dose_edep.raw", dtype=np.float64) 

# Reshape en 3D pour correspondre à DimSize
data = data.reshape((dim_x, dim_y, dim_z))

# Afficher la matrice 8x8 de dose déposée
print("Dose déposée dans le microarray :")
print(data[:, :, 0])  # On affiche uniquement la 1ère couche (z=0)

plt.imshow(data[:, :, 0], cmap='inferno', origin='lower')
plt.colorbar(label="Dépôt d'énergie ")
plt.title("Distribution de la dose dans le microarray")
plt.xlabel("X (voxel)")
plt.ylabel("Y (voxel)")
plt.show()




#  Charger un fichier ROOT
filename = "test_hits.root"

#  Ouvrir le fichier ROOT
with uproot.open(filename) as file:
    #  Lister les objets disponibles dans le fichier
    print("Objets disponibles dans le fichier ROOT :")
    obj = file["Hit_collec"]
    print(f"Type de 'Hit_collec' : {obj.classname}")
    print(f"Branches disponibles dans {"Hit_collec"} :")
    print(obj.keys())
    for key in file.keys():
        print(f" - {key}")

    branch_name = "TotalEnergyDeposit"  # nom de la branche à regarder
    data = obj[branch_name].array(library="np")
    
    
    #  Affichage
    import matplotlib.pyplot as plt
    plt.hist(data, bins=50, color='blue', alpha=0.7, edgecolor='black')
    plt.xlabel(branch_name)
    plt.ylabel("Occurrences")
    plt.title(f"Distribution de {branch_name}")
    plt.show()



# units
"""m = gate.g4_units.m
cm = gate.g4_units.cm
mm = gate.g4_units.mm
keV = gate.g4_units.keV
MeV = gate.g4_units.MeV
Bq = gate.g4_units.Bq

#volume
volume = sim.add_volume("BoxVolume","firstvolume")
volume.size = [4*cm, 4*cm, 4*cm]
volume.material = "G4_WATER"

#source X-ray  -
incident_beam = sim.add_source('GenericSource' , 'MyFirstSource')
incident_beam.position.type = "point"
incident_beam.direction.type = "iso"
incident_beam.n = 10
incident_beam.energy.mono = 69 * keV

stats = sim.add_actor("SimulationStatisticsActor","stats")
stats.track_types_flag = True
stats.steps_per_event_flag = True

print("lancement de la simulation")
sim.run()
print("fin de la simulation")

print(stats)
     
print("nombre total de particules primaires émises :", stats.counts.events)


fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection= '3d')

#afficher sur un histogramme le nombre de steps par tracks
tab = np.zeros((2,10)) #chaque colonne correspond à une track et la deuxième ligne contient le nombre de steps d'une track donnée
for event_index , steps_per_event in stats.counts.steps_per_event.item() :
    tab.append([event_index, steps_per_event])


print(tab) """   