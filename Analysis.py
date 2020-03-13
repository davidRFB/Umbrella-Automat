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

colvar = list(np.float_(open("COLVAR.dat").read().splitlines()))
frames = list(np.float_(open("frames.dat").read().splitlines()))

CV_i = colvar[-1]
CV_f = colvar[0]
#0,529177 convert bohr to amnstrong15
C_ts = (CV_i -CV_f)*0.5291773/2
dE_tot = C_ts*2


E_ts_guess = float(input("\n guess value for E_ts \n"))

K_Force  = 6*5*E_ts_guess/(C_ts)**2
Tem = float(input("\n Temperature Value \n"))
beta = 1/((const.Boltzmann*const.Avogadro/4184)*(Tem))
# %%
dE_wind= 2/math.sqrt(K_Force*beta)  


## calculate the number of the windows for the umbrella sampling
Num_Wind = math.ceil(2*C_ts/dE_wind)  

## iterate over directories 

directories=[d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and "CV" in d ]

cv_vals = []
os.system("mkdir Analysis")

## 
directories.sort(key=lambda x: float(x[2:]))


    

for j,i in enumerate(directories):
    try:
        CV_dat1 = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(i,i))
        CV_dat2 = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(directories[j+1],directories[j+1]))
        #histogram_intersection(CV_data[:,1],CV_dat2[:,1])
        plt.figure()
        plt.title("{}".format("{}/{}-COLVAR.metadynLog".format(i,i)))
        y_1,x_1, _ =plt.hist(CV_dat1[:,1],bins=30)
        
        min_a1 = min(x_1)
        max_a1 = max(x_1)
        middle_a1 = x_1[round(len(x_1)/2)]
        plt.show()
        plt.figure()
        plt.title("{}/{}-COLVAR.metadynLog".format(directories[j+1],directories[j+1]))
        y_2,x_2, _ =plt.hist(CV_dat2[:,1],bins=30)
        #print(y_2.max(),x_2.max(),x_1[round(len(x_1)/2)])
        min_b = min(x_2)
        max_b = max(x_2)
        middle_b = x_2[round(len(x_2)/2)]
        plt.show()
        #print(y_2.max(),x_2.max(),x_1[round(len(x_1)/2)])
        print("resolution {}".format (2*(middle_a1-middle_b)/((max_a1-min_a1)-(max_b-min_b))))
    except:
        print("out")

##################plots generation 
# %%
## create concatenated histogram 
plt.figure(figsize=(15,10))
for i in directories:
    CV_data = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(i,i))
    plt.title(" Histogram of all CVs  ".format(i))
    plt.hist(CV_data[:,1],bins=10,alpha=0.5)
    plt.axvline(np.mean(CV_data[:,1]),linestyle="--")

plt.savefig("./Analysis/Full_histrograms.pdf")    
plt.show()

## generate all Cvs vs frames 
plt.figure(figsize=(10,10))
for i in directories:
    CV_data = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(i,i))
    plt.title("Colvar {} vs time ".format(i))
    plt.plot(CV_data[:,0],CV_data[:,1])
    plt.axhline(np.mean(CV_data[:,1]),linestyle="--")
plt.savefig("./Analysis/all_Cvs_time.pdf")    
plt.show()
# %%
data_new = []
os.system("mkdir W")
for k,i in enumerate(directories):
    CV_data = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(i,i))[:,1]
    K_sinf = np.genfromtxt("./test_ks.txt")
    with open("w{}".format(k),"+a") as wind_file:
        for j in range(len(CV_data)):
            data_new.append((CV_data[j],K_sinf[k][0],K_sinf[k][1]))       
            wind_file.write("{} {} {} \n".format(CV_data[j],K_sinf[k][0],K_sinf[k][1]))


# %%

### generate every CV plot 
'''
for i in dir_CV:
    plt.figure()
    CV_data = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(i,i))
    plt.title("Colvar {} vs time ".format(i))
    plt.plot(CV_data[:,0],CV_data[:,1])
    plt.savefig("./Analysis/Colvar{}_vs_time ".format(i))
    plt.figure()
    plt.title("Colvar {} Histogram  ".format(i))
    plt.hist(CV_data[:,1],bins=50)
    plt.savefig("./Analysis/Colvar{}_Histogram".format(i))
'''
# %%

def histogram_intersection(h1, h2):
    sm = 0
    for i in range(len(h1)):
        sm += min(h1[i], h2[i])
    return sm
'''
for j,i in enumerate(dir_CV):
    try:
        CV_data = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(i,i))
        print("{}/{}-COLVAR.metadynLog".format(dir_CV[j+1],dir_CV[j+1]))
        CV_dat2 = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(dir_CV[j+1],dir_CV[j+1]))
        histogram_intersection(CV_data[:,1],CV_dat2[:,1])
    except:
        print("out")
#print(y,x[(len(x)-1)],x)
#print(np.bincount(CV_data).argmax())
#y,x = np.histogram(CV_data[:,1],bins=21)  
#print(y,len(x)/2)
#print(y[int(len(y)/3)]/2,x[int(len(x)/2)])
#closest_val = min(enumerate(colvar), key=lambda x: abs(x[1]-i))
#plt.figure()
#plt.hist(CV_data[:,1],bins=20)
# %%
CV_dat1 = np.genfromtxt("075CV/075CV-COLVAR.metadynLog")
CV_dat2 = np.genfromtxt("060CV/060CV-COLVAR.metadynLog")
    

y_1,x_1 = np.histogram(CV_dat1[:,1],bins=21)
y_2,x_2 = np.histogram(CV_dat2[:,1],bins=21)
    
plt.hist(CV_dat1[:,1],bins=20)  
plt.hist(CV_dat2[:,1],bins=20)      
    
    
    

def histogram_intersection(h1, h2):
    sm = 0
    for i in range(13):
        sm += min(h1[i], h2[i])
    return sm

print(histogram_intersection(CV_dat1[:,1],CV_dat2[:,1]))    
    
    
#closest_val = min(enumerate(y_2), key=lambda x: abs(x[1]-histogram_intersection(CV_dat1[:,1],CV_dat2[:,1])))    
'''    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




