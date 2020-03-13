#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 18:13:23 2020

@author: david
"""
# %%
import numpy as np  
import sys, os
import math 
import scipy.constants as const
import mdtraj as mdt
import time 

## read file from Colective variable from the Steered Molecular dynamic (SMD)
Colvar_data_file = sys.argv[1]

# create COVLAR and frames files 
os.system("cat "+Colvar_data_file+" | awk '{ print $2}' > COLVAR.dat")
os.system("cat "+Colvar_data_file+" | awk '{ print $1}' > frames.dat")

colvar = list(np.float_(open("COLVAR.dat").read().splitlines()))
frames = list(np.float_(open("frames.dat").read().splitlines()))

CV_i = colvar[-1]
CV_f = colvar[0]
#0,529177 convert bohr to amnstrong
C_ts = (CV_i -CV_f)*0.5291773/2
dE_tot = C_ts*2


E_ts_guess = float(input("\n guess value for E_ts \n"))

K_Force  = 6*5*E_ts_guess/(C_ts)**2
Tem = float(input("\n Temperature Value \n"))
beta = 1/((const.Boltzmann*const.Avogadro/4184)*(Tem))

dE_wind= 2/math.sqrt(K_Force*beta)


## calculate the number of the windows for the umbrella sampling
Num_Wind = math.ceil(2*C_ts/dE_wind)  


print(" Reading and proccesing trayectory with mdtraj ------- \n\n")
## load the trayectory in mdtraj to extract frames
start_time = time.time()
#traj = mdt.load_xyz("magl_new_steer-pos-1.xyz",top="./MAGL_ARA.prmtop")
print("--- %s seconds of processing time \n " % (time.time() - start_time))
## Script writing

print(" Writing files ------- \n\n")
K_s = []
for i in np.linspace(CV_f,CV_i,Num_Wind):
    cv = i
    # creates directory of each window
    try:
        os.system("mkdir CV{:.3f}".format(i))
    except:
        print("directorio creado")
    new_direc_cv="CV{:.3f}".format(i)
    #base_dir = os.getcwd()
    #os.chdir(new_direc_cv)
    closest_val = min(enumerate(colvar), key=lambda x: abs(x[1]-i))
    traj[int(round(frames[closest_val[0]]/10))].save_pdb("{}/{:.3f}_CV.pdb".format(new_direc_cv,i))
    if(abs(i) > 0.6*(C_ts)):
        K_Force = 4/(beta*(dE_wind + 0.4*dE_wind)**2)
    else:
        K_Force= 6*5*E_ts_guess/(C_ts)**2
    #print(i,K_Force)
    K_s.append(K_Force)
    with open("test_.txt","a+") as test_ks:
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
      parm_file_name glicerol.top
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
  MM_INDEX   1824
  MM_INDEX   1826
        MM_INDEX   1827
  MM_INDEX   4134
  MM_INDEX   4137
  MM_INDEX   4140
  MM_INDEX   4143
  MM_INDEX   3654
  MM_INDEX   3657
    &END QM_KIND
    &QM_KIND N
  MM_INDEX   4138
  MM_INDEX   4142
   &END QM_KIND
    &QM_KIND O
  MM_INDEX   1819
  MM_INDEX   1825
  MM_INDEX   3658
  MM_INDEX   3659
  MM_INDEX   115453
    &END QM_KIND
    &QM_KIND H
  MM_INDEX   1822
  MM_INDEX   1823
  MM_INDEX   1828
  MM_INDEX   1829
  MM_INDEX   1831
  MM_INDEX   1832
  MM_INDEX   4135
  MM_INDEX   4136
  MM_INDEX   4139
  MM_INDEX   4141
  MM_INDEX   4144
  MM_INDEX   3655
  MM_INDEX   3656
  MM_INDEX   115454
  MM_INDEX   115455
    &END QM_KIND

   &LINK
      ALPHA 1.40
         LINK_TYPE IMOMM
      MM_INDEX 1815
      QM_INDEX 1818
      RADIUS 0.80
    &END LINK
    &LINK
      ALPHA 1.40
      LINK_TYPE IMOMM
      MM_INDEX 1830
      QM_INDEX 1827
      RADIUS 0.80
    &END LINK
        &LINK
      ALPHA 1.40
      LINK_TYPE IMOMM
      MM_INDEX 3652
      QM_INDEX 3654
      RADIUS 0.80
    &END LINK
    &LINK
      ALPHA 1.40
      LINK_TYPE IMOMM
      MM_INDEX 4132
      QM_INDEX 4134
      RADIUS 0.80
    &END LINK

  &END QMMM

  &SUBSYS
    &CELL
      ABC 144.825   93.890  130.731
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      CONN_FILE_NAME glicerol.top
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
        ATOMS 115453 1824 1824 1819
        COEFFICIENT -1.0
      &END
    &END

  &END SUBSYS
&END FORCE_EVAL

&GLOBAL
  PRINT_LEVEL LOW
  PROJECT nombre
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
        umb_cv.writelines("        K  "+str("%.4f" % (K_Force*0.000446253514585234))+"\n")
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

