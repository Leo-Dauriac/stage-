#!/usr/bin/env python
import argparse
import os
import json
from multiprocessing import Pool

import matplotlib.pyplot as plt
import opengate as gate
import uproot
import numpy as np
import pyvista as pv
from scipy.optimize import curve_fit

from scipy.spatial.transform import Rotation as R
import spekpy as sp
import critical_angle 
#import psutil

DEFAULT_OUTPUT = 'output' #file name where single root file are stocked
DEFAULT_NUMBER_OF_JOBS = int(1e1)
DEFAULT_NUMBER_OF_PARTICLES = int(1e8)
DEFAULT_POSITION = [0, 0, 0]
DEFAULT_STEP_SIZE = float(1)
DEFAULT_WORLD = "G4_Galactic" #default material of the world volume
DEFAULT_SAMPLE = "G4_WATER"   #default material of the sample volume
DEFAULT_FILE_NAME = "simu_test"


def opengate_run(
    output=DEFAULT_OUTPUT,
    job_id=0,
    number_of_particles=DEFAULT_NUMBER_OF_PARTICLES,
    visu=True, 
    verbose=False,
    polarization=False,
    energy_type=False,
    world=DEFAULT_WORLD,
    sample_material=DEFAULT_SAMPLE,
    File_name=DEFAULT_FILE_NAME,
    position=DEFAULT_POSITION,
     ):

    if verbose:
        print(f"Running MC… " + str(locals())) #locals = variable locales soit les paramètres de la fonction

    if polarization:
        polarization = [1, 0, 0]
        #print("Polarization is enabled")
        #print(f"{polarization}")
    else:
        polarization = [0, 0, 0]  
        #print("Polarization is disabled") 
        #print(f"{polarization}")



    output_base = output #output_format(air_length_mm, water_length_mm, output)
    #output_base est un chemin
    new_pos = position
    

    # Units
    sec = gate.g4_units.second
    m = gate.g4_units.m
    cm = gate.g4_units.cm
    mm = gate.g4_units.mm
    um = gate.g4_units.um
    keV = gate.g4_units.keV
    
    # Simulation
    
    sim = gate.Simulation()

    sim.random_engine = 'MersenneTwister'
    sim.random_seed = 'auto'
    #sim.run_timing_intervals = [[0 * sec, 1 * sec]]
    sim.check_volumes_overlap = False
    sim.visu = visu 
    sim.visu_type = 'vrml'
    sim.g4_verbose = False
    sim.progress_bar = verbose
    sim.number_of_threads = 1

    # Misc

    yellow = [1, 1, 0, 1]
    blue = [0, 0, 1, 1]
    red = [1, 0, 0, 1]

    # Geometry
    #sim.volume_manager.add_material_database(gate.utility.get_tests_folder() / '..' / 'data' / 'GateMaterials.db')
    sim.volume_manager.add_material_database("/home/ldauriac/Software/opengate/opengate/opengate/tests/data/GateMaterials.db")
    sim.world.material = world
    print(f"World's material is {sim.world.material}")
    #sim.world.material = 'Air'
    #sim.world.material = "G4_He"
    #sim.world.material = "G4_Galactic"
    sim.world.size = [40 * cm, 40 * cm, 40 * cm]
    sim.world.color = [0, 0, 0, 0]

    # Phantom

    #Polyethylene film volume where the sample is placed
    membrane = sim.add_volume("BoxVolume" , "membrane")
    membrane.size = [26*mm, 76*mm , 1.4*um]
    membrane.translation = [0, 0, 0]
    membrane.material = "G4_POLYETHYLENE" 
    rotation_matrix_obj = R.from_euler('y', -15, degrees=True).as_matrix()
    membrane.rotation = rotation_matrix_obj

    cell_vol = sim.add_volume("TubsVolume" , "cell_vol") #Cell volume 
    cell_vol.rmin = 0
    cell_vol.rmax = 10*mm
    cell_vol.dz = 20*um
    cell_vol.material = "G4_WATER"
    cell_vol.translation = [0,0,22*um]
    rotation_matrix_cell = R.from_euler('zy', [90,-15], degrees=True).as_matrix()
    cell_vol.rotation = rotation_matrix_cell
    cell_vol.color = red

    sample = sim.add_volume("TubsVolume" , "sample") #target volume simulating the fluorescence emitter concentration
    sample.rmin = 0
    sample.rmax = 100*um
    sample.dz = 20*um
    sample.translation = [0,0,22*um]
    sample.material = sample_material #material of the sample, can be changed with the command line argument
    print(f"Sample's material is {sample.material}")
    #sample.material = "Water_3mgI"
    #sample.material = "H2O_I_Gd_mix"
    #sample.material = "Water_50mgI" #50mg/g d'Iode, on peut potentiellement monter jusqu'à 300mg/g comme du ioméron 400mg/ml de densité 1.35
    #sample.material = "Water_300mgI"
    #sample.material = "ppmPtry"
    rotation_matrix_sample = R.from_euler('zy', [90,15], degrees=True).as_matrix()
    sample.rotation = rotation_matrix_sample
    sample.color = blue

    """collimator = sim.add_volume("TubsVolume" , "collimator")
    collimator.rmin = 5*mm
    collimator.rmax = 7*mm
    collimator.dz = 1*cm
    collimator.material = "G4_Pb"
    collimator.translation = [-5.8*cm, 0, 0]
    rotation_matrix_detector = R.from_euler('y', -90, degrees=True).as_matrix()
    collimator.rotation = rotation_matrix_detector 
    collimator.color = red

    pin_hole = sim.add_volume("TubsVolume", "pin_hole")
    pin_hole.rmin = 3*mm
    pin_hole.rmax = 5*mm
    pin_hole.dz = 1*mm
    pin_hole.material = "G4_Pb"
    pin_hole.translation = [-4.9*cm, 0, 0]
    rotation_matrix_detector = R.from_euler('y', -90, degrees=True).as_matrix()
    pin_hole.rotation = rotation_matrix_detector 
    pin_hole.color = blue"""
    

    # Beam
    ebin = 0.5
    energy_poly = 120
    spek = sp.Spek(kvp=energy_poly, th=12, dk=ebin,physics="casim", mu_data_source="pene", mas=1.0)
    spek.filter('Ag', 1)
    (x,y) = spek.get_spectrum()
    xb = x - ebin/2.
    kV=[energy_poly]
    xb = np.append(xb,kV)

    source = sim.add_source('GenericSource' , "mybeam")
    source.particle = "gamma"
    source.position.type = "point"
    source.direction.type = "momentum"
    source.direction.momentum = [0, 0, -1]
    source.polarization = polarization
    #source.polarization = [1, 0, 0]
    source.position.translation = new_pos #décalée sur le pixel le plus à gauche pour commencer (elle prendra +100 sur x à chaque trans)
    source.n = number_of_particles

    if energy_type:
        print("Beam is in polychromatic mode")
        source.energy.type = "spectrum_histogram" 
        source.energy.spectrum_energy_bin_edges = xb * keV
        source.energy.spectrum_weights = y
    else: 
        print("Beam is in monochromatic mode")
        source.energy.mono = 100 * keV    
    
    
    # Physics list

    sim.physics_manager.physics_list_name = "G4EmLivermorePolarizedPhysics"
    sim.physics_manager.em_parameters.fluo = True
    sim.physics_manager.em_parameters.auger = True  
    sim.physics_manager.em_parameters.auger_cascade = True
    sim.physics_manager.global_production_cuts.gamma = 0.001 * um
    sim.physics_manager.global_production_cuts.electron = 100 * um

    # HPGe Detector + Hits Collection
    
    HPGe = sim.add_volume("TubsVolume" , "HPGe")
    HPGe.rmin = 0
    HPGe.rmax = 5*mm
    HPGe.dz = 3*mm
    HPGe.material = "G4_Ge"
    HPGe.color = yellow
    HPGe.translation = [-6.5*cm, 0, 0] 
    rotation_matrix_detector = R.from_euler('y', 90, degrees=True).as_matrix()
    HPGe.rotation = rotation_matrix_detector 

    

    #actor
    HPGe_detection = sim.add_actor("DigitizerHitsCollectionActor" , "HPGe_test_detection")
    HPGe_detection.attached_to = ["HPGe"]
    HPGe_detection.attributes = ["TotalEnergyDeposit"]
    #HPGe_detection.output_filename = os.path.join(output_base, f"HPGE_detection_{job_id}.root")
    HPGe_detection.output_filename = os.path.join(output_base, f"{File_name}_{job_id}.root")
    #HPGe_detection.output_filename = os.path.join(output_base, f"Rasterscan_test_{job_id}.root")
    parfilter = sim.add_filter("ParticleFilter", "parfilter")
    parfilter.particle = "gamma"
    HPGe_detection.filters.append(parfilter)

    # Particle stats

    stat = sim.add_actor('SimulationStatisticsActor', 'stat')
    stat.track_types_flag = True
    stat.output_filename = os.path.join(output_base, f'stats_{job_id}.txt')

    sim.run()
    
    #visualization
    # === Charger le fichier VRML avec PyVista ===
    """print("📂 Chargement du fichier VRML pour affichage...")
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
    pl.show()  # Afficher la scène 3D"""

def root_visu(
    output=DEFAULT_OUTPUT,
    verbose=False,
    File_name=DEFAULT_FILE_NAME,
):
    def print_verbose(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)
    try:
        #data = uproot.concatenate(os.path.join(output, 'HPGE_detection_*.root'), library='np')
        data = uproot.concatenate(os.path.join(output, f'{File_name}_*.root'), library='np')
        #va chercher tous les fichiers qui commencent par HPGE_detection_ dans le dossier output et rassemble les données en une structure de type NumPy 
        #data = uproot.concatenate(os.path.join(output, 'Rasterscan_test_*.root'), library='np')
        ws = 1000 * data['TotalEnergyDeposit'] #multiplie tte les données par 1000 pour passer de MeV à keV
        
    except FileNotFoundError:
        print("file not found")

    #myfile = uproot.recreate("HPGe.root")
    myfile = uproot.recreate(f"{File_name}.root")
    #myfile["HPGe"] = {"TotalEnergyDeposit": ws} #met dans le Ttree "HPGE" une branche TotalEnergyDeposit contenant les données de ws
    myfile[f"{File_name}"] = {"TotalEnergyDeposit": ws}
    #myfile = uproot.recreate("Rasterscan.root")
    #myfile["Rasterscan"] = {"TotalEnergyDeposit": ws}

    #showing final concatenated root file    
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel("Photon energy [keV]")
    ax.set_ylabel("Counts")
    
    plt.hist(ws, bins=500, range=(0,100))

    plt.savefig(os.path.join(output, 'histo.pdf'))
    plt.show()

    delta1 = critical_angle.critical_angle(100)[0]
    delta2 = critical_angle.critical_angle(100)[1]
    theta = critical_angle.critical_angle(100)[2]
    print(f"delta1 : {delta1}")
    print(f"delta2 : {delta2}")
    print(f"critical angle : {theta}")

def opengate_pool_run(
    output,
    number_of_jobs, #mettre autant que de positions sinon on ne simule pas toutes les particules
    number_of_particles,
    visu,
    verbose,
    polarization,
    energy_type,
    world,
    sample_material,
    File_name,
    #step,
    #size,
):

    with Pool(maxtasksperchild=1) as pool:

        results = []
        
        
        number_of_particles_per_job = int(number_of_particles / number_of_jobs)
        mm = gate.g4_units.mm
        um = gate.g4_units.um
        cm = gate.g4_units.cm
        position = [0, 0, 20*cm] #position initiale à 30um -> trop proche de la source
        #delta_x = 0.1*mm #pas du raster scan pour passer d'un pixel à l'autre
        #delta_x = step*mm
        #delta_y = -step*mm
        job_id = 0
        #image_size = size

        """ for y in range(image_size) : #lignes
            for x in range(image_size) : #colonnes
                job_id += 1
                 # Calcule la nouvelle position à partir de la base
                position = [
                    x * delta_x,
                    y * delta_y,
                    position[2]
                ] 
                #le carré décrit est balayé de gauche à droite et de haut en bas
                
                copied_position = position.copy()
                print(f"launching job #{job_id}/{number_of_jobs} with position x={copied_position[0]/mm:.2f} mm, y={copied_position[1]/mm:.2f} mm")

                result = pool.apply_async(opengate_run, kwds={
                'output': output,
                'job_id': job_id, 
                'number_of_particles': number_of_particles_per_job,
                'visu': visu,
                'verbose': verbose,
                'position': copied_position
                #'position' : position 
                })
                results.append(result)"""



        
        for i in range(number_of_jobs):
            
                job_id = i+1
                #copied_position = position.copy()
                print(f"launching job #{job_id}/{number_of_jobs}")
                #print(f"launching job #{job_id}/{number_of_jobs} with position x={copied_position[0]/mm:.2f} mm")
                result = pool.apply_async(opengate_run, kwds={
                'output': output,
                'job_id': job_id, 
                'number_of_particles': number_of_particles_per_job,
                'visu': visu,
                'verbose': verbose,
                'position' : position,
                'polarization': polarization,
                'energy_type': energy_type,
                'world' : world,
                'sample_material' : sample_material,
                'File_name': File_name
                #'position' : copied_position,
                })
                results.append(result)
                #position[0] += delta_x

            
        

        pool.close()
        pool.join()

        for result in results:
            result.wait()
            if not result.successful():
                print("Failure in MC simulation")
                exit(1)




    root_visu(
        output=output,
        verbose=verbose,
        File_name=File_name,
    )

                
def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--output', help="Path of outputs", default=DEFAULT_OUTPUT)
    parser.add_argument('-j', '--number-of-jobs', help="Number of jobs", default=DEFAULT_NUMBER_OF_JOBS, type=int)
    parser.add_argument('-n', '--number-of-particles', help="Number of generated particles (total)", default=DEFAULT_NUMBER_OF_PARTICLES, type=int)
    parser.add_argument('--visu', help="Visualize Monte Carlo simulation", default=False, action='store_true')
    parser.add_argument('--verbose', '-v', help="Verbose execution", default=False, action='store_true')
    parser.add_argument('-p', '--polarization', help="Polarization of the beam", default=False, action='store_true')
    parser.add_argument('--energy_type', help="Type of energy distribution", default=False, action='store_true')#True active le mode polychromatique, rendre la commande plus clair après
    parser.add_argument('-w', '--world', help="World's material", default=DEFAULT_WORLD, type=str)#pas besoin de mettre les guillemets au moment de la commande
    parser.add_argument('-s', '--sample_material', help="sample's material", default=DEFAULT_SAMPLE, type=str)
    parser.add_argument('-f', '--File_name', help="name of the output file", default=DEFAULT_FILE_NAME, type=str)
    #parser.add_argument('-step', '--step', help="size of a step between two succesive positions en mm", default=DEFAULT_STEP_SIZE, type=float)
    #parser.add_argument('-size', '--size', help="number of pixels in lines/columns", default=DEFAULT_STEP_SIZE, type=int)
    #parser.add_argument('-pos', '--number-of-positions', help="Number of positions to simulate", default=DEFAULT_NUMBER_OF_POSITIONS, type=int)
    #argument position initiale ?     
    #argument pour le nombre de position à couvrir (taille de l'image à balayer en gros)
    args_info = parser.parse_args()

    opengate_pool_run(**vars(args_info))
    #opengate_run(output, 1, number_of_particles, visu, verbose)
if __name__ == '__main__':
    main()