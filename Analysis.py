#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:28:40 2020

@author: david
"""
import numpy as np  
import sys, os
import math 
import scipy.constants as const
import matplotlib.pyplot as plt

print("######  This script automatize the process of the umbrella sampling analysis  ###### \n ")
print("###### from the SMD_to_US.py script. Generating histograms plots and corrected US   ###### \n")
print("######           should be excuted as $ python Analysis.py                    ###### ")



colvar = list(np.float_(open("COLVAR.dat").read().splitlines()))
frames = list(np.float_(open("frames.dat").read().splitlines()))

CV_i = colvar[-1]
CV_f = colvar[0]

#0,529177 used to convert bohr to amnstrong

C_ts = (CV_i -CV_f)*0.5291773/2
dE_tot = C_ts*2
print("Values of E_ts_guess must be the same of the SMD_to_US execution \n")
E_ts_guess = float(input("\n guess value for E_ts (Kj/mol) \n"))

K_Force  = 6*5*E_ts_guess/(C_ts)**2
Tem = float(input("\n Temperature Value  (K) \n"))
beta = 1/((const.Boltzmann*const.Avogadro/4184)*(Tem))
# %%
# Window width by equation (18) in 
# Kästner, J., & Thiel, W. (2006). Analysis of the statistical error in umbrella sampling simulations by umbrella integration.
#  The Journal of Chemical Physics, 124(23), 234106. doi:10.1063/1.2206775 
dE_wind= 2/math.sqrt(K_Force*beta)  


## calculate the number of the windows for the umbrella sampling
Num_Wind = math.ceil(2*C_ts/dE_wind)  

## iterate over directories 
directories=[d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and "CV" in d ]
cv_vals = []
os.system("mkdir Analysis")

## organize directories in order of CV (e.i [min_Cv...Max_Cv])
directories.sort(key=lambda x: float(x[2:]))

for j,i in enumerate(directories):
    try:
        # CV data
        CV_dat1 = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(i,i))
        CV_dat2 = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(directories[j+1],directories[j+1]))
        #histogram of CV1
        y_1,x_1 =np.histogram(CV_dat1[:,1],bins=bins_num)
        #parameters of CV1 histogram
        min_a1 = min(x_1)
        max_a1 = max(x_1)
        middle_a1 = x_1[round(len(x_1)/2)]
        #histogram of CV2
        y_2,x_2 =np.histogram(CV_dat2[:,1],bins=bins_num)
        #parameters of CV2 histogram
        min_b = min(x_2)
        max_b = max(x_2)
        middle_b = x_2[round(len(x_2)/2)]
        #Resolution definition (Chromatography definition)
        resultion_2_1 = 2*(middle_b-middle_a1)/((max_b-min_b)+(max_a1-min_a1))
        print("resolution {} between {} and {} \n".format (resultion_2_1,i,directories[j+1] ))
        if(resultion_2_1> 0.85):
            print("check values for CVs {} and {} \n".format(i,directories[j+1]))
    except:
        print("out")

##################plots generation 
# %%
## create concatenated histogram 
plt.figure(figsize=(15,10))
for i in directories:
    CV_data = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(i,i))
    plt.title(" Histogram of all CVs  ")
    plt.hist(CV_data[:,1],bins=10,alpha=0.5)
    plt.axvline(np.mean(CV_data[:,1]),linestyle="--")
    plt.legend(loc=(1.05,0.15))
plt.savefig("./Analysis/Full_histrograms.jpg")    
plt.show()

## generate all Cvs vs frames 
plt.figure(figsize=(10,10))
for i in directories:
    CV_data = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(i,i))
    plt.title("Colvar {} vs time ".format(i))
    plt.plot(CV_data[:,0],CV_data[:,1])
    plt.axhline(np.mean(CV_data[:,1]),linestyle="--")
    plt.legend(loc=(1.05,0.05))
plt.savefig("./Analysis/all_Cvs_time.jpg")    
plt.show()

### generate every CV plot 

for i in directories:
    plt.figure()
    CV_data = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(i,i))
    plt.title("Colvar {} vs time ".format(i))
    plt.plot(CV_data[:,0],CV_data[:,1])
    plt.savefig("./Analysis/Colvar{}_vs_time ".format(i))
    plt.figure()
    plt.title("Colvar {} Histogram  ".format(i))
    plt.hist(CV_data[:,1],bins=50)
    plt.savefig("./Analysis/Colvar{}_Histogram".format(i))

    
###### PMF analysis ###### 
data_new = []
os.system("mkdir W")
for k,i in enumerate(directories):
    CV_data = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(i,i))[:,1]
    K_sinf = np.genfromtxt("./test_ks.txt")
    with open("w{}".format(k),"+a") as wind_file:
        for j in range(len(CV_data)):
            data_new.append((CV_data[j],K_sinf[k][0],K_sinf[k][1]))       
            wind_file.write("{} {} {} \n".format(CV_data[j],K_sinf[k][0],K_sinf[k][1]))



    
    
    
    
    
    
    
    
    
    
    
    




