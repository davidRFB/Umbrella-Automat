#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 18:13:23 2020

@author: David Ricardo Figueroa Blanco dr.figueroa10@uniandes.edu.co
"""
# %%

# modules !!
import numpy as np  
import sys, os
import math 
import scipy.constants as const
import mdtraj as mdt
import time 

print("###### This script automatize the process of the umbrella sampling execution  ###### \n ")
print("######  from a Steered Molecular Dynamic (SMD) previously performed in CP2K   ###### \n")
print("######  should be excuted as $ python SMD_to_US.py [COLVAR-FILE] [.xyz trayectory]   ###### ")

## read file from Colective variable from the Steered Molecular dynamic (SMD)
Colvar_data_file = sys.argv[1]

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
traj = mdt.load_xyz("magl_new_steer-pos-1.xyz",top="./MAGL_ARA.prmtop")
print("--- %s seconds of processing time \n " % (time.time() - start_time))
## Script writing

print(" ----- Writing files ------- \n\n")
K_s = []
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

    traj[int(round(frames[closest_val[0]]/10))].save_pdb("{}/{:.3f}_CV.pdb".format(new_direc_cv,i))
    
    ## Definition of a narrow K constant far from the CV of the Ts 
    # Actual definition ( 60 % of the defined constant)
    if(abs(i) > 0.6*(C_ts)):
        K_Force = 4/(beta*(dE_wind + 0.4*dE_wind)**2)
    else:
        K_Force= 6*5*E_ts_guess/(C_ts)**2
    #print(i,K_Force)
    K_s.append(K_Force)
    ## creation of a file that will be used in the creation of a Potential Mean Force (PMF) profile
    with open("pmf_data.txt","a+") as test_ks:
        test_ks.write("{:.4f} {:.4f} \n".format(K_Force*0.000446253514585234,i))
        
# %%
    start_time = time.time()
    with open ('{}/US_{:.3f}CV.sh'.format(new_direc_cv,i),"w") as umb_cv:
        umb_cv.write('''&FORCE_EVAL
  METHOD QMMM

  &DFT
   CHARGE -1
   BASIS_SET_FILE_NAME BASIS_MOLOPT
   POTENTIAL_FILE_NAME GTH_POTENTIALS
!   RESTART_FILE_NAME magl_ptos_steer-RESTART.wfn

    &MGRID
      CUTOFF 320
      NGRIDS 5
      COMMENSURATE
    &END MGRID

    &QS
      EPS_DEFAULT 1.0E-12
    &END QS

    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 5.0E-7
      MAX_SCF 100
      &OUTER_SCF
        EPS_SCF 5.0E-7
        MAX_SCF 10
      &END
      &OT
        MINIMIZER DIIS
!       PRECONDITIONER FULL_SINGLE_INVERSE
!       PRECONDITIONER FULL_KINETIC
        PRECONDITIONER FULL_ALL
        ENERGY_GAP 0.001
      &END
    &END SCF

    &XC
      &XC_FUNCTIONAL BLYP
      &END XC_FUNCTIONAL
      &vdW_POTENTIAL
        DISPERSION_FUNCTIONAL PAIR_POTENTIAL
        &PAIR_POTENTIAL
           TYPE DFTD3
           CALCULATE_C9_TERM .FALSE.
           REFERENCE_C9_TERM .FALSE.
           LONG_RANGE_CORRECTION .FALSE.
           PARAMETER_FILE_NAME dftd3.dat
!          VERBOSE_OUTPUT
           REFERENCE_FUNCTIONAL BLYP
           EPS_CN 0.01
           R_CUTOFF 20.0
        &END PAIR_POTENTIAL
      &END vdW_POTENTIAL
      &XC_GRID
        XC_SMOOTH_RHO NN10
        XC_DERIV SPLINE2_SMOOTH
      &END XC_GRID
    &END XC

    &POISSON
      POISSON_SOLVER WAVELET
      PERIODIC NONE
      &WAVELET
        SCF_TYPE 60
      &END
    &END

  &PRINT
   &E_DENSITY_CUBE
    &EACH
     MD 500
    &END
    ADD_LAST NUMERIC
    STRIDE 2 2 2
   &END
  &END

  &END DFT
    &MM
    &FORCEFIELD
      &SPLINE
       EMAX_SPLINE 10.0
       RCUT_NB 10.
      &END
      parm_file_name ../MAGL_ARA.prmtop 
      parmtype AMBER
      VDW_SCALE14 0.5
      EI_SCALE14  0.8333333
    &END FORCEFIELD
    &POISSON
      &EWALD
        EWALD_TYPE spme
        ALPHA .24
        GMAX 145  94  131
        RCUT 10.
        O_SPLINE 6
      &END EWALD
    &END POISSON
  &END MM

  &QMMM
    USE_GEEP_LIB 10
    &CELL
      ABC 21. 21. 21.
      PERIODIC NONE
    &END CELL
    E_COUPL GAUSS
    &INTERPOLATOR
      EPS_R 1.0e-14
      EPS_X 1.0e-14
      MAXITER 200
    &END INTERPOLATOR

 &QM_KIND C
  MM_INDEX   1818
  MM_INDEX   3603
  MM_INDEX   3606
  MM_INDEX   4083
  MM_INDEX   4086
  MM_INDEX   4089
  MM_INDEX   4092
  MM_INDEX   75600
  MM_INDEX   75601
  MM_INDEX   75602
  MM_INDEX   75605
  MM_INDEX   75606
  MM_INDEX   75607
    &END QM_KIND
    &QM_KIND N
  MM_INDEX   4087
  MM_INDEX   4091
    &END QM_KIND
    &QM_KIND O
  MM_INDEX   1821
  MM_INDEX   3607
  MM_INDEX   3608
  MM_INDEX   75603
  MM_INDEX   75608
  MM_INDEX   75609
  MM_INDEX   75610
    &END QM_KIND
    &QM_KIND H
  MM_INDEX   1819
  MM_INDEX   1820
  MM_INDEX   1822
  MM_INDEX   3604
  MM_INDEX   3605
  MM_INDEX   4084
  MM_INDEX   4085
  MM_INDEX   4088
  MM_INDEX   4090
  MM_INDEX   4093
  MM_INDEX   75629
  MM_INDEX   75630
  MM_INDEX   75631
  MM_INDEX   75632
  MM_INDEX   75635
  MM_INDEX   75636
  MM_INDEX   75637
  MM_INDEX   75638
  MM_INDEX   75639
  MM_INDEX   75640
  MM_INDEX   75641
    &END QM_KIND

  &LINK
      ALPHA 1.40
      LINK_TYPE IMOMM
      MM_INDEX 1816
      QM_INDEX 1818
      RADIUS 0.80
    &END LINK
    &LINK
      ALPHA 1.40
      LINK_TYPE IMOMM
      MM_INDEX 3601
      QM_INDEX 3603
      RADIUS 0.80
    &END LINK
    &LINK
      ALPHA 1.40
      LINK_TYPE IMOMM
      MM_INDEX 4081
      QM_INDEX 4083
      RADIUS 0.80
    &END LINK
    &LINK
      ALPHA 1.40
      LINK_TYPE IMOMM
      MM_INDEX 75594
      QM_INDEX 75600
      RADIUS 0.80
    &END LINK


  &END QMMM

  &SUBSYS
    &CELL
      ABC 144.825   93.890  130.731
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      CONN_FILE_NAME ../MAGL_ARA.prmtop
      CONNECTIVITY AMBER \n''')
        umb_cv.writelines("      COORD_FILE_NAME   {:.3f}_CV.pdb \n".format(i))
        umb_cv.write('''
      COORDINATE PDB
      PARA_RES .FALSE. 
      &DUMP_PDB
      &END
    &END TOPOLOGY

    &KIND C
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-BLYP-q4
    &END KIND
    &KIND N
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-BLYP-q5
    &END KIND
    &KIND O
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-BLYP-q6
    &END KIND
    &KIND H
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-BLYP-q1
    &END KIND

    &COLVAR
      &DISTANCE_FUNCTION
        ATOMS 75603 75602 75602 1821
        COEFFICIENT -1.0
      &END
    &END

  &END SUBSYS
&END FORCE_EVAL

&GLOBAL
  PRINT_LEVEL LOW \n ''')
        umb_cv.writelines("      PROJECT  CV{:.3f} \n".format(i))
        umb_cv.write('''
  RUN_TYPE MD
  PREFERRED_FFT_LIBRARY FFTSG
&END GLOBAL
&MOTION

  &MD
   ENSEMBLE NVT
   STEPS 3000
   TIMESTEP 0.5
   TEMPERATURE 300.0
   &THERMOSTAT
    TYPE CSVR
     REGION DEFINED
     &DEFINE_REGION
       MM_SUBSYS ATOMIC
     &END
     &DEFINE_REGION
       QM_SUBSYS ATOMIC
     &END
    &CSVR
      TIMECON 50
     &END
    &END
   &END MD

  &PRINT
   &TRAJECTORY
    FORMAT XMOL
    &EACH
     MD 10
    &END
   &END TRAJECTORY
   &RESTART
    &EACH
     MD 100
    &END
    ADD_LAST NUMERIC
   &END RESTART
   &FORCES OFF
    &EACH
     MD 100
    &END
      ADD_LAST NUMERIC
   &END
   &VELOCITIES OFF
   &END
  &END PRINT

   &FREE_ENERGY
     &METADYN
       LANGEVIN  F
       DO_HILLS  F
       LAGRANGE  F
       NT_HILLS  1
       WW     0.
       TEMPERATURE 1000
       &METAVAR
         SCALE   1.
         COLVAR  1
       &END METAVAR
       &PRINT
         &COLVAR
           COMMON_ITERATION_LEVELS  3
           &EACH
             MD  10
           &END EACH
         &END COLVAR
       &END PRINT
     &END METADYN
   &END FREE_ENERGY

  &CONSTRAINT
    &COLLECTIVE
      COLVAR 1
      INTERMOLECULAR \n '''
      )
        umb_cv.writelines("        TARGET   "+str("%.4f" % i)+"\n")
        umb_cv.write('''!      TARGET_GROWTH 0.0001
      &RESTRAINT
      ''')
        umb_cv.writelines("   K  "+str("%.4f" % (K_Force*0.000446253514585234))+"\n")
        umb_cv.write('''      &END
     &END COLLECTIVE
  &END

&END MOTION

!&EXT_RESTART
!  RESTART_FILE_NAME magl_ptos_steer-1.restart
!  RESTART_THERMOSTAT .FALSE.
&END
'''
                )
        print("--- %s seconds writing  \n " % (time.time() - start_time))
#files = os.listdir()               
#for i in files:
#    if (".sh" in i ):
#        print(i)

