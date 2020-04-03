#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:28:40 2020

@author: David Ricardo Figueroa Blanco 
@email:dr.figueroa10@uniandes.edu.co
"""

import numpy as np  
import sys, os
import math 
import scipy.constants as const
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pathlib2 import Path


print("######  This script automatize the process of the umbrella sampling analysis  ###### \n ")
print("###### from the SMD_to_US.py script. Generating histograms plots and corrected US   ###### \n")
print("######           should be excuted as $ python Analysis.py                    ###### ")

message = None

try:
    CV_i = float(os.popen("grep Cv_i Report_SMD_to_US.txt | awk ' {print $3}'").read())
    CV_f = float(os.popen("grep Cv_f Report_SMD_to_US.txt | awk ' {print $3}'").read())
    E_ts_guess =  float(os.popen("grep E_gs Report_SMD_to_US.txt | awk ' {print $3}'").read())
    Tem= float(os.popen("grep Temp Report_SMD_to_US.txt | awk ' {print $3}'").read())
except :
    print("Report_SMD_to_US.txt not found \n ")

#0,529177 used to convert bohr to amnstrong

C_ts = (CV_i -CV_f)*0.5291773/2
dE_tot = C_ts*2

K_Force  = 6*5*E_ts_guess/(C_ts)**2

beta = 1/((const.Boltzmann*const.Avogadro/4184)*(Tem))

# Window width by equation (18) in 
# Kästner, J., & Thiel, W. (2006). Analysis of the statistical error in umbrella sampling simulations by umbrella integration.
#  The Journal of Chemical Physics, 124(23), 234106. doi:10.1063/1.2206775 
dE_wind= 2/math.sqrt(K_Force*beta)  

Ks_data = np.genfromtxt("./pmf_data.txt")


## calculate the number of the windows for the umbrella sampling
Num_Wind = math.ceil(2*C_ts/dE_wind)  

## iterate over directories 
directories=[d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and "CV" in d ]
cv_vals = []
os.system("mkdir Analysis")

## organize directories in order of CV (e.i [min_Cv...Max_Cv])
directories.sort(key=lambda x: float(x[2:]))
bins_num=60
resolution =[]
for j,i in enumerate(directories):
    #try:
        # CV data
        if(j==len(directories)-1):
            print("Final")
            break
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
        resolution.append(resultion_2_1)
        
        with open ('Report_Analysis_US.txt',"a+") as report:
            report.write("resolution {} between {} and {} \n".format (resultion_2_1,i,directories[j+1] ))
            
        if(resultion_2_1> 1):
            with open ('Report_Analysis_US.txt',"a+") as report:
                report.write("checking values for CVs {} and {} \n".format(i,directories[j+1]))
            print("checking values for CVs {} and {} \n".format(i,directories[j+1]),j)
            #check if previous CV is nor overlaping 
            if(resolution[-2]>1.0):
                tochangeCV=directories[resolution.index(resolution[-1])]
                try:
                    path = Path("{}/{}US.sh".format(i,i))
                    text = path.read_text()
                    text = text.replace("K  "+str(Ks_data[:,0][j]), "K  "+str(Ks_data[:,0][j]*0.5))
                    text = text.replace("!   RESTART_FILE_NAME toedit_byresolution", "   RESTART_FILE_NAME {}-RESTART.wfn".format(i))
                    text = text.replace("PROJECT  {}".format(i),"PROJECT  {}_2k".format(i))
                    path.write_text(text)
                    with open ('Report_Analysis_US.txt',"a+") as report:
                         report.write("changin K of "+str(Ks_data[:,0][j])+" by "+str(Ks_data[:,0][j]*0.5)+"in {} \n".format(tochangeCV))
                except:
                    print("nofile")

##################plots generation 
# %%
## create concatenated histogram 
plt.figure(figsize=(15,10))
for i in directories:
    CV_data = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(i,i))
    plt.title(" Histogram of all CVs  ")
    plt.hist(CV_data[:,1],bins=10,alpha=0.5,label="{}".format(i))
    plt.axvline(np.mean(CV_data[:,1]),linestyle="--")
    plt.legend(loc=(1.05,0.15))
plt.savefig("./Analysis/Full_histrograms.jpg",bbox_inches="tight")    
#plt.show()

## generate all Cvs vs frames 
plt.figure(figsize=(10,10))
for i in directories:
    CV_data = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(i,i))
    plt.title("Colvar {} vs time ".format(i))
    plt.plot(CV_data[:,0],CV_data[:,1],label="{}".format(i))
    plt.axhline(np.mean(CV_data[:,1]),linestyle="--")
    plt.legend(loc=(1.05,0.05))
plt.savefig("./Analysis/all_Cvs_time.jpg",bbox_inches="tight")    
#plt.show()

### generate every CV plot 

for i in directories:
    plt.figure()
    CV_data = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(i,i))
    plt.title(" {} vs time ".format(i))
    plt.plot(CV_data[:,0],CV_data[:,1],label="{}".format(i))
    plt.savefig("./Analysis/{} vs time .png".format(i[:-1]))

    plt.figure()
    plt.title(" {} Histogram  ".format(i))
    plt.hist(CV_data[:,1],bins=50)
    plt.savefig("./Analysis/{} histogram.png".format(i[:-1]))

#%%
###### PMF analysis ###### 
data_new = []
try:
    os.system("mkdir W")
except:
    print("W directory already created !!! ")
## create files w with Colvar to make the pmf analysis.
for k,i in zip(range(len(directories)-1),directories):
    CV_data = np.genfromtxt("{}/{}-COLVAR.metadynLog".format(i,i))[:,1]
    K_sinf = np.genfromtxt("./pmf_data.txt")
    with open("./W/w{}".format(k),"+a") as wind_file:
        for j in range(len(CV_data)):
            data_new.append((CV_data[j],K_sinf[k][0],K_sinf[k][1]))       
            wind_file.write("{:.6f} {:.6f} {:.6f} \n".format(CV_data[j],K_sinf[k][0],K_sinf[k][1]))

os.system("mkdir pmf_test")
#os.system("mv W pmf_test")
os.system("cp job  pmf_test ")
os.system("cp umbrella_integration.x pmf_test" )

print("Execute job inside pmf_test")

#%%

pmf_final = np.genfromtxt("./pmf_test/fe_ui.xy",skip_header=2)
plt.figure()
#converting to kcal/mol
pmf_final_kcal = (pmf_final[:,1])*627.5
plt.title("Energy Profile")
plt.xlabel(" Reaction coordinate ")
plt.ylabel(" Energy (kcal/mol) ")
plt.plot(pmf_final[:,0],pmf_final_kcal,label=" $E_a = $ {:.4f}".format(max(pmf_final_kcal)))
plt.legend(fontsize=16,loc="best")
plt.savefig("pmf.pdf")
plt.show()

# %%
