import opengate as gate
from opengate.tests import utility
from scipy.spatial.transform import Rotation
import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import math
import xraylib


z_tag_con = np.zeros((16,2))
z_tag_con = [ ("Y89", 1.25), 
        ("In113", 1.25),
        ("La139", 2.5), 
        ("Ce140", 1.25), 
        ("Pr141", 5), 
        ("Nd150", 1.25), 
        ("Sm147", 0.625), 
        ("Eu153", 0.15), 
        ("Gd158", 0.4), 
        ("Tb159", 2.5), 
        ("Dy164", 5), 
        ("Ho165", 0.6), 
        ("Er170", 0.32), 
        ("Tm169", 1.2), 
        ("Yb173", 0.125), 
        ("Lu175", 1.25), ] 

def choix_element() :
    e = np.random.uniform(0,1)
    if e <= 0.051 :
     return "G4_Y"
    elif e <= 0.102 : 
     return "G4_In"
    elif e <= 0.205 : 
     return "G4_La"
    elif e <= 0.256 :
     return "G4_Ce"
    elif e<= 0.456 : 
     return "G4_Pr"
    elif e <= 0.507 :
     return "G4_Nd"
    elif e <= 0.0532 :
     return "G4_Sm"
    elif e <= 0.539 :
     return "G4_Eu"
    elif e <= 0.555 :
     return "G4_Gd"
    elif e <= 0.657 : 
     return "G4_Tb"
    elif e <= 0.857 : 
     return "G4_Dy"
    elif e <= 0.881 :
     return "G4_Ho"
    elif e <= 0.894 : 
     return "G4_Er"
    elif e <= 0.944 :
     return "G4_Tm"
    elif e <= 0.949 :
     return "G4_Yb"
    else : 
     return "G4_Lu"
    

   # Fonction pour vérifier s'il y a chevauchement
def is_overlapping(new_pos, existing_positions, min_distance):
    for pos in existing_positions:
        distance = np.linalg.norm(np.array(new_pos) - np.array(pos))
        if distance < min_distance:
            return True  # Chevauchement détecté
    return False  # Pas de chevauchement




import numpy as np

def random_point_in_cylinder(parent_radius, parent_height, child_radius):
    """
    Génère une position aléatoire pour un z_tag dans un sample cylindrique, en évitant l'overlap.
    - parent_radius : rayon du sample (ex: 0.5 mm)
    - parent_height : hauteur du sample (ex: 10 µm = 0.01 mm)
    - child_radius : rayon du z_tag (ex: 0.05 mm)
    """

    max_radius = parent_radius - child_radius  # Évite que z_tag touche les bords
    max_height = parent_height - 2 * child_radius  # Évite que z_tag sorte du cylindre

    # Coordonnées cylindriques : rayon et angle
    r = max_radius * np.sqrt(np.random.uniform(0, 1))  # Distribution uniforme en r²
    theta = np.random.uniform(0, 2 * np.pi)  # Angle entre 0 et 2π

    # Coordonnées cartésiennes dans le disque (base du cylindre)
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    # Position aléatoire en hauteur
    z = np.random.uniform(-max_height / 2, max_height / 2)  # Centré autour de 0

    return [x, y, z]


#def random_point_in_sphere(parent_radius, child_radius):
    """
    Génère une position aléatoire pour un z_tag dans un sample, en évitant l'overlap.
    - parent_radius : rayon du sample (ex: 0.5 mm)
    - child_radius : rayon du z_tag (ex: 0.05 mm)
    """
    max_radius = parent_radius - child_radius  # Évite que z_tag touche les bords
    r = max_radius * np.cbrt(np.random.uniform(0, 1))
    theta = np.arccos(1 - 2 * np.random.uniform(0, 1))  # [0, π]
    phi = np.random.uniform(0, 2 * np.pi)  # [0, 2π]

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return [x, y, z]




def DD_prob(nbr_parti_init) :
  Z_Pr = 59
  Z_Si = 14
  
  energy = 69 #keV
  shell = xraylib.K_SHELL #couche concernée
  L1_K = xraylib.KL1_LINE #transition concernée
  fluo_energy = xraylib.LineEnergy(Z_Pr, L1_K) #énergie de transition
  
  µrho_Pr = xraylib.CS_Photo(Z_Pr, energy) #g.cm⁻2
  rho_Pr = xraylib.ElementDensity(Z_Pr)
  µ_Pr = µrho_Pr * rho_Pr #coefficient d'atténuation linéique par effet PE dans le Pr
 
  µrho_Si = xraylib.CS_Photo(Z_Si, fluo_energy) #g.cm⁻2
  rho_Si = xraylib.ElementDensity(Z_Si)
  µ_Si = µrho_Si * rho_Si #coefficient d'atténuation linéique par effet PE dans le Si
  
  
  microarray_prob = 1-(np.exp(-µ_Pr*1e-4)) #probabilité d'intéraction par effet PE dans le microarray
  Pr_fluo_yield = xraylib.FluorYield(Z_Pr, shell) #rendement de fluorescence
  solid_angle_prob = (((2*(np.pi*(3.1**2)))/30**2))/(4*np.pi) #probabilité d'angle solide
  detector_prob = 1-(np.exp(-µ_Si*5*1e-4)) #probabilité d'intéraction par effet PE dans le détecteur
  prob_tot = microarray_prob*Pr_fluo_yield*solid_angle_prob*detector_prob #probabilité total d'intéraction avec le détecteur
  detected_particle = prob_tot*nbr_parti_init 
  
  return [prob_tot, detected_particle]





