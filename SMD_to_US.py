#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 18:13:23 2020

@author: David Ricardo Figueroa Blanco dr.figueroa10@uniandes.edu.co
"""
import numpy as np  
import sys, os
import math 
import scipy.constants as const
import mdtraj as mdt
import matplotlib.pyplot as plt
from pathlib2 import Path
import argparse
import time
import shutil


# Instantiate the parser
parser = argparse.ArgumentParser(description='Automatation of Umbrella Sampling Calculation')

parser.add_argument("-c",dest='SMD_COLVAR',  required=True,
                    help='Colective variable data <.metadyLog > ')
parser.add_argument("-x",dest='traj_xyz', required=True ,
                    help='SMD trayectory from cp2k <.xyz > ')
parser.add_argument("-p",dest='topology_file', required=True, 
                    help='topology file amber <.top .mrtop> ')
parser.add_argument("-t",dest='template_cp2k', required=True, 
                    help='cp2k input <Template US > ')


args = parser.parse_args()

print(args.SMD_COLVAR)
print(args.traj_xyz)
print(args.topology_file)

print("###### This script automatize the process of the umbrella sampling execution  ###### \n ")
print("######  from a Steered Molecular Dynamic (SMD) previously performed in CP2K   ###### \n")
print("######  should be excuted as $ python SMD_to_US.py [COLVAR-FILE] [.xyz trayectory]   ###### ")

## read file from Colective variable from the Steered Molecular dynamic (SMD)
Colvar_data_file = args.SMD_COLVAR
SMD_traj = args.traj_xyz
topol_file = args.topology_file
template = args.template_cp2k
# create COVLAR and frames files 
os.system("cat "+Colvar_data_file+" | awk '{ print $2}' > COLVAR.dat")
os.system("cat "+Colvar_data_file+" | awk '{ print $1}' > frames.dat")

colvar = list(np.float_(open("COLVAR.dat").read().splitlines()))
frames = list(np.float_(open("frames.dat").read().splitlines()))


CV_i = colvar[-1]
CV_f = colvar[0]

#0,529177 used to convert bohr to amnstrong
C_ts = (CV_i -CV_f)*0.5291773/2
dE_tot = C_ts*2


E_ts_guess = float(input("\n guess value for E_ts (kJ/mol) \n"))

K_Force  = 6*5*E_ts_guess/(C_ts)**2
Tem = float(input("\n Temperature Value (K) \n"))
## calculation of beta factor for window determination
beta = 1/((const.Boltzmann*const.Avogadro/4184)*(Tem))


# Window width by equation (18) in 
# Kästner, J., & Thiel, W. (2006). Analysis of the statistical error in umbrella sampling simulations by umbrella integration.
#  The Journal of Chemical Physics, 124(23), 234106. doi:10.1063/1.2206775 
dE_wind= 2/math.sqrt(K_Force*beta)


## calculate the number of the windows for the umbrella sampling
Num_Wind = math.ceil(2*C_ts/dE_wind)  


print(" ----- Reading and proccesing trayectory with mdtraj ------- \n\n")
## load the trayectory in mdtraj to extract frames
start_time = time.time()
#traj = mdt.load_xyz(SMD_traj,top=topol_file)
print("--- %s seconds of processing time \n " % (time.time() - start_time))
## Script writing

print(" ----- Writing files ------- \n\n")
K_s = []
try:
    os.remove("pmf_data.txt")
except:
    print("first execution ! Good luck ! \n ")


for i in np.linspace(CV_f,CV_i,Num_Wind):
    # creates directory of each window
    try:
        os.system("mkdir CV{:.3f}".format(i))
    except:
        print("directorio creado")
    new_direc_cv="CV{:.3f}".format(i)

    #find the closest value of Collective variable (CV) into the COLVAR file from the SMD
    closest_val = min(enumerate(colvar), key=lambda x: abs(x[1]-i))

    ## extract pdb from the trajectory

    #traj[int(round(frames[closest_val[0]]/10))].save_pdb("{}/CV{:.3f}.pdb".format(new_direc_cv,i))
    
    ## Definition of a narrow K constant far from the CV of the Ts 
    # Actual definition ( 60 % of the defined constant)
    if(abs(i) > 0.6*(C_ts)):
        K_Force = 4/(beta*(dE_wind + 0.4*dE_wind)**2)
    else:
        K_Force= 6*5*E_ts_guess/(C_ts)**2
    #print(i,K_Force)
    K_s.append(K_Force)
    ## creation of a file that will be used in the creation of a Potential Mean Force (PMF) profile
    # Handle errors while calling os.remove()

    with open("pmf_data.txt","a+") as test_ks:
        test_ks.write("{:.4f} {:.4f} \n".format(K_Force*0.000446253514585234,i))
        
    
    shutil.copy(template,"{}/CV{:.3f}US".format(new_direc_cv,i))
    
    
    path = Path("{}/CV{:.3f}US".format(new_direc_cv,i))
    text = path.read_text()
    text = text.replace("topoledit", "../{}".format(topol_file))
    text = text.replace("coordfiledit", "CV{:.3f}.pdb".format(i))
    text= text.replace("nametoedit","CV{:.3f}".format(i))
    text= text.replace("cvvaltoedit","{}".format(i))
    text= text.replace("constantktoedit","{:.4f}".format(K_Force*0.000446253514585234))
    path.write_text(text)
    #write job to execution 
    with open("CV{:.3f}/job".format(i),"w") as job:
        job.write('''#!/bin/bash
#SBATCH --time=160:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH --partition=medium
#SBATCH --err=myJob.err 
#SBATCH --out=myJob.out 
#SBATCH --mail-type=ALL \n''')
        job.writelines("#SBATCH --job-name=CV{:.3f}".format(i))
        job.write('''module load cp2k/4.1-ssmp \n''')
        job.writelines("cp2k.ssmp -i CV{:.3f} -o output \n".format(i)) 


try:
    os.remove("Report_SMD_to_US.txt")
except:
    print("first execution ! Good luck ! \n ")


with open ('Report_SMD_to_US.txt',"w") as report:
        report.write(time.ctime())
        report.write("\n")
        report.write(''' This report contains information related with the Umbrella sampling calculation ''')
        report.write("\n")
        report.writelines("Cv_i =   {} \n".format(CV_i))
        report.writelines("Cv_f =   {} \n".format(CV_f))
        report.writelines("Temp =   {} \n".format(Tem))
        report.writelines("E_gs =   {} \n".format(E_ts_guess))
        report.writelines("N_win =   {} \n".format(Num_Wind))
        report.writelines("Wind_Dist= {:.4f} \n".format(dE_wind))
        report.writelines("K = {:.4f} \n".format(K_Force))
        report.writelines("K_2= {:.4f} \n".format(6*5*E_ts_guess/(C_ts)**2))
        report.write("QM Region \n")

os.system("sed -n '/QM_KIND C/,/LINK/p' {}/CV{:.3f}US | head -n -2 >> Report_SMD_to_US.txt ".format(new_direc_cv,i))


try:
    os.remove("dftd3.dat")
    os.remove("BASIS_MOLOPT")
    os.remove("GTH_POTENTIALS")

except:
    print("first execution ! Good luck ! \n ")

print("Downloading the Basis set and GTH potential files from https://github.com/misteliy/cp2k/tree/master/tests/QS  ")
os.system("wget https://raw.githubusercontent.com/misteliy/cp2k/master/tests/QS/dftd3.dat ")
os.system("wget https://raw.githubusercontent.com/misteliy/cp2k/master/tests/QS/BASIS_MOLOPT")
os.system("wget https://raw.githubusercontent.com/misteliy/cp2k/master/tests/QS/GTH_POTENTIALS")

os.system("for i in *CV*/; do cd $i ; echo sbatch job ; cd .. ; done")