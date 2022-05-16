#!/usr/bin/env python -u

import subprocess

import numpy as np
import os


import argparse
parser = argparse.ArgumentParser(description='Process some model coupling arguments.')
parser.add_argument('-t','--dt', dest='timestep', type=float, default=1.0,
                    help='Define the coupling time step')
parser.add_argument('-ys','--ys', dest='startyear', type=float, default=-123.0,
                    help='Define the startyear')
parser.add_argument('-ye','--ye', dest='endyear', type=float, default=0.0,
                    help='Define the endyear')
parser.add_argument('-it','--it', dest='iteration', type=int, default=1,
                    help='Number of iteration')
parser.add_argument('-p','--cpupertask', dest='cpupertask', type=int, default=1,
                    help='Define the number of CPU per task')
parser.add_argument('-m', '--module', dest="module", default='prepare',
                    choices=['prepare', 'pism2vilma', 'vilma2pism'],
                    help='Name of coupler modue to call')



args = parser.parse_args()
print(args)
print("\nRun coupling."+args.module+"() with "+str(args.timestep)+" kyr timestep with "+str(args.cpupertask)+" CPU")


### tool that passes command strings to bash ########
def run(cmd):
  print(cmd)
  return subprocess.check_call(cmd, shell=True)



### set start, end and time interval #################
secperyear=365.0*24.0*3600.0

# number of glacial cycle runs with updated topography
iteration=args.iteration
timestep=args.timestep
cpupertask=args.cpupertask

# define time intervals in kyr
inityear=-246.0
btime = args.startyear
etime = btime + timestep
ftime = args.endyear

# pism-related restart and regrid files ############################
pismoutlast='pism/results/pism_last.nc'
pismout4vilma='pism/results/pism4vilma.nc'
pismoutfile='pism/results/paleo.nc'
pisminfile='pism/results/paleo_inp.nc'

# PISM readable set of initial conditions at ys
if iteration==1:
  pismstart='data/boot_16km_tw.nc'
else:
  pismstart='data/plig_16km_it'+str(iteration-1)+'.nc'

# make ice5/6G file cdo readable ################################
icegcdo="data/ice6g/ice6g_cdo_246.nc"
timestepyr=str(int(timestep*1000))+"years"

icegdt="data/ice6g/ice6g_"+timestepyr+"_2gc.nc"

# combination of ICE6g and Bedmap2 Antarctica (topg0.nc) on Gauss Legendre gris n128
#ncks -A -v topo ../02_vilma_standalone/data/ice6g/bedmap2ice6g217.nc data/topo0.nc
topoorig = "data/topo"+str(iteration-1)+".nc"
topgpismorig = "pism/results/topg"+str(iteration-1)+".nc"

# VILMA gloabal grid size and time steps, ice and ocean densities, consistent with PISM
load_pism_file = "data/ice6g/loadh.inp"

# defines how structure is read in: resolution, load history, reference topography and initial ice distribution
load_hist_file = "inp_ice6g_pismant/load_hist.inp"
load_hist_init="data/ice6g/Ice0.nc"

#########################################################################################
########################################################################################



def prepare():

  # as in vilma run script
  if os.path.exists("io.tmp"):
    os.remove("io.tmp")

  # pism output data
  if not os.path.exists("pism/results"):
    os.makedirs("pism/results")

  # vilma output data
  if not os.path.exists("out"):
    os.makedirs("out")
  os.system("touch out/vega_oce.dat out/vega_rpt.dat out/vega_deg1.dat")

  if not os.path.exists("out/topo0.nc"):
    os.system("ncks -v topo "+topoorig+" out/topo0.nc")

  # vilma restart files
  dir_restart = "restart"
  if not os.path.exists(dir_restart):
    os.system('mkdir '+dir_restart)

  # make a copy of the 3d Earth structure fields
  #os.system('cp data/visc3d_new.nc out/visc3d_new.nc')
  os.system('ln -s ../data/visc3d_new.nc out/visc3d_new.nc')

  # initfile for first PISM iteration
  if not os.path.exists(pisminfile):
    run("cp "+pismstart+" "+pisminfile)
    intime=str(float(int(inityear*1000.0*secperyear)))
    run("ncap2 -O -s 'time(:)={"+intime+"}' "+pisminfile+" "+pisminfile)

  # first snapshot in PISM histroy passed to VILMA
  if not os.path.exists(pismoutlast):
    run("ncks -A -v thk,mask "+pismstart+" "+pismoutlast)
    intime=str(float(int(inityear*1000.0*secperyear)))
    run("ncap2 -O -s 'time(:)={"+intime+"}' "+pismoutlast+" "+pismoutlast)

  # prepare original ice load data to have equidistant time steps
  if not os.path.exists(icegdt):
    print("\nInterpolate for Ice5/6G..")
    run("cdo -P "+str(cpupertask)+" inttime,-246000-01-01,00:00:00,"+timestepyr+" "+icegcdo+" "+icegdt)

  # initial bed topography in PISM
  if not os.path.exists(topgpismorig):
    run("ncks -A -v topg "+pismstart+" "+topgpismorig)
    run("ncrename -v topg,topg0 "+topgpismorig)
    ##run("ncks -O -s 'time=-3878928000000' "+topgpismorig+" "+topgpismorig)

  print("\nRun PISM.. ########################################################################")



def pismtovilma():

        btimew=str(np.around(btime,decimals=3))
        etimew=str(np.around(etime,decimals=3))

        # convert to years
        dtyr=int(timestep*1000)
        btimeyr=str(int(np.around(btime*1000/dtyr,decimals=1)*dtyr))
        etimeyr=str(int(np.around(etime*1000/dtyr,decimals=1)*dtyr))

        # convert to seconds
        btimes=str(float(int(btime*1000.0*secperyear)))
        etimes=str(float(int(etime*1000.0*secperyear)))

        vilma2pismgridlast="out/vilma2pism"+btimew.zfill(6)+".nc"
        if np.float(btime)!=inityear and not os.path.exists(vilma2pismgridlast):
          print("\n"+vilma2pismgridlast+" does not exits, exit iterations...")
          exit(1)

        
        #### remap pism thk to vilma grid ############################################################
        # merge last and new PISM snapshot to define history (more intermediate snaphots needed?)
        print("\nPrepare PISM data for VILMA forcing..")
        #os.system("cdo mergetime -selname,thk,mask "+pismoutlast+" pism/results/paleo.nc "+pismoutfile)
        if os.path.exists(pismout4vilma): 
          os.remove(pismout4vilma)
        run("ncks -A -v thk,mask "+pismoutfile+" "+pismout4vilma)
        run("cdo -P "+str(cpupertask)+" -O mergetime "+pismoutlast+" "+pismout4vilma+" "+pismout4vilma)

        # delete PISM history
        run("ncatted -a history_of_appended_files,global,d,c, \
                           -a history,global,d,c, "+pismout4vilma)

        # rename pism output files
        run("rm "+pismoutlast)
        run("ncks -A -v thk,mask "+pismoutfile+" "+pismoutlast)
        
        # select two snapshots from ice5g history and merge pism output
        print("\nPrepare Ice5/6G data for VILMA forcing..")
        iceg_input="out/ice6g_"+timestepyr+etimew.zfill(6)+".nc"
        run("cdo -P "+str(cpupertask)+" selyear,"+btimeyr+","+etimeyr+" "+icegdt+" "+iceg_input)

        # remap from PISM (stere) grid to VILMA (gaussian) grid
        print("\nRemap PISM to Ice5/6G data..")
        pism2vilmagrid="pism/results/pism2vilma"+etimew.zfill(6)+".nc"
        # FIXME: in order to save time, save remap weights!
        #run("cdo -P "+str(cpupertask)+" remapbil,"+ice5g_input+" "+pismout4vilma+" "+pism2vilmagrid)
        run("cdo -P "+str(cpupertask)+" remapbic,"+iceg_input+" "+pismout4vilma+" "+pism2vilmagrid)

        # merge PISM output variables to ice5g history
        run("ncatted -a _FillValue,mask,d,, -a missing_value,mask,d,, "+pism2vilmagrid)
        run("ncks -A -v thk,mask "+pism2vilmagrid+" "+iceg_input)

        # merge PISM output thickness into Antarctic ice5g history
        run("ncap2 -O -s 'where(mask>0) Ice=thk' "+iceg_input+" "+iceg_input)

        # delete PISM output variables from ice5g history
        iceg_temp=iceg_input+"-tmp"
        run("ncks -O -x -v thk,mask "+iceg_input+" "+iceg_temp)
        run("mv "+iceg_temp+" "+iceg_input)

        # make time variable VILMA readable
        run("ncap2 -O -s 'time(:)={"+btimew+","+etimew+"}' "+iceg_input+" "+iceg_input)
        run("ncrename -v time,epoch -d time,epoch "+iceg_input)
        run("ncatted -a units,epoch,o,c,'ka BP' \
                     -a long_name,epoch,o,c,'Epoch' \
                     -a calendar,epoch,o,c,'proleptic_gregorian' "+iceg_input)

        # save initial ice thickness as reference for all VILMA runs
        if btime==inityear:
            run("ncks -d epoch,0 -v Ice "+load_hist_init+" out/ice6g_ref.nc")
            run("ncap2 -O -s 'epoch(:)={"+str(btime)+"}' out/ice6g_ref.nc out/ice6g_ref.nc")
            run("ncks -O --mk_rec_dmn epoch out/ice6g_ref.nc out/ice6g_ref.nc")
            run("ncks -d epoch,1 "+iceg_input+" out/ice6g_step1.nc")
            run("rm "+iceg_input)
            run("ncrcat out/ice6g_ref.nc out/ice6g_step1.nc "+iceg_input)
            run("rm out/ice6g_step1.nc")

        # copy ice5g_input to data and edit inp_pism/load_hist.inp
        load_hist_txt = load_pism_file+"\n"+iceg_input+"\n"+topoorig+" "+load_hist_init
        with open(load_hist_file, "w") as f:
          f.write(load_hist_txt)        
 
        # update start and end year for VILMA run in inp_pism/wepochs.inp
        load_pism_txt = "256 512 910 1028\n2\n"+btimew+"\n"+etimew
        #load_pism_txt = "512 1024 910 1028\n2\n"+btimew+"\n"+etimew
        with open(load_pism_file, "w") as f:
          f.write(load_pism_txt)

        ### run vilma #############################################################################
        print("\nRun VILMA.. ########################################################################")


def vilmatopism():


        btimew=str(np.around(btime,decimals=3))
        etimew=str(np.around(etime,decimals=3))

        # convert to years
        dtyr=int(timestep*1000)
        btimeyr=str(int(np.around(btime*1000/dtyr,decimals=1)*dtyr))
        etimeyr=str(int(np.around(etime*1000/dtyr,decimals=1)*dtyr))

        # convert to seconds
        btimes=str(float(int(btime*1000.0*secperyear)))
        etimes=str(float(int(etime*1000.0*secperyear)))

        vilma2pismgrid="out/vilma2pism"+etimew.zfill(6)+".nc"
        run("cdo -P "+str(cpupertask)+" remapbil,"+pismout4vilma+" out/rsl.nc "+vilma2pismgrid)
        #FIXME: save remap weights to speed up coupling

        # get only last time slice and save
        try:
          run("ncks -A -v rsl -d time,1 "+vilma2pismgrid+" "+vilma2pismgrid+"-tmp")
          run("mv "+vilma2pismgrid+"-tmp "+vilma2pismgrid)
        except:
          #run("cdo tinfo "+vilma2pismgrid)
          print(vilma2pismgrid+" has only one time slice")


        run("ncks -A -v topg0 "+topgpismorig+" "+vilma2pismgrid)
        run("ncap2 -O -s 'topg=topg0-rsl' "+vilma2pismgrid+" "+vilma2pismgrid)

        # prepare PISM restart (with modified topg) with next coupling time step
        run("cp "+pismoutfile+" "+pisminfile)
        run("ncks -A -v topg "+vilma2pismgrid+" "+pisminfile)
        run("ncap2 -O -s 'time(:)={"+etimes+"}' "+pisminfile+" "+pisminfile)
          
        # add rsl to pism file for diagnostic
        run("ncks -A -v rsl "+vilma2pismgrid+" "+pismoutfile)
        
        # add dbdt=(rsl-rsl_previous)/dtyr for diagnostic
        try:
          run("ncdiff -v rsl "+pismoutfile+" pism/results/paleo"+btimew.zfill(6)+".nc out/dbdt.nc")
          #cdo sub 
        except:
          run("ncks -O -v rsl "+pismoutfile+" out/dbdt.nc") #in case of -123.0 ka BP
        run("ncap2 -O -s 'dbdt=-1000.0*rsl/"+str(dtyr)+"' out/dbdt.nc out/dbdt.nc")        
        run("ncatted -O -a units,dbdt,o,c,'mm year-1' \
                        -a long_name,dbdt,o,c,'bedrock uplift rate' \
                        -a standard_name,dbdt,o,c,'tendency_of_bedrock_altitude' out/dbdt.nc")
        run("ncks -A -v dbdt out/dbdt.nc "+pismoutfile)
        run("rm out/dbdt.nc")
        run("ncap2 -O -s 'time(:)={"+etimes+"}' "+pismoutfile+" "+pismoutfile)

        # make a copy of pism output file for the record, avery 100 years
        run("ncatted -O -a history,global,d,, "+pismoutfile)
        run("mv "+pismoutfile+" pism/results/paleo"+etimew.zfill(6)+".nc")

        run("ncatted -O -a history,global,d,, pism/results/ts_paleo.nc")
        run("mv pism/results/ts_paleo.nc pism/results/ts_paleo"+etimew.zfill(6)+".nc")

        # make a copy of vilms rsl output file for the record
        run("mv out/rsl.nc out/rsl"+etimew.zfill(6)+".nc")

        #end of iteration with coupling time step #########################################

        print("######################################################################################\n\n")



if args.module == 'prepare':
    prepare()
elif args.module == 'pism2vilma':
    pismtovilma()
elif args.module == 'vilma2pism':
    vilmatopism()
else:
    print("Something went wrong in the choice of coupler modules "+args.module)

