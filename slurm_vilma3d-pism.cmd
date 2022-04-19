#!/bin/bash

#SBATCH --qos=medium

#SBATCH --job-name=vilmapism_def3d_dt100_it1
#SBATCH --account=pism
#SBATCH --output=./log-vilmapism-%j.out
#SBATCH --error=./log-vilmapism-%j.err

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=albrecht@pik-potsdam.de

#SBATCH --time=6-22:59:00

#SBATCH --constraint=haswell,tasksmax
#SBATCH --exclusive
#SBATCH --no-requeue


### PISM #################################################

#SBATCH --ntasks=64
#SBATCH --tasks-per-node=16

#SBATCH packjob


### VILMA3d and python coupler ###########################

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

#export OMP_NUM_THREADS=16
#export KMP_AFFINITY=verbose,granularity=thread,compact,1

#export SLURM_CPUS_PER_TASK=$OMP_NUM_THREADS
export CPUS_PER_TASK=16

#echo coupler:
#echo $SLURM_NTASKS
#echo $SSLURM_CPUS_PER_TASK

# Some attributes are inherited (e.g. the second invocation gets the same 5 minute timelimit)
# https://slurm.schedmd.com/heterogeneous_jobs.html
###########################################################################################

# specify iteration
it=1

# specify coupling time steps
timeint=.1
ys=-246.0
ye=0

# run coupler ##########################################################################################################
# prepare coupler
srun --pack-group=1 python run_coupler.py -m "prepare" --dt=${timeint} -ys ${ys} -ye ${ye} -it ${it} --cpupertask=${CPUS_PER_TASK}

export LC_NUMERIC="en_US.UTF-8"
pad=1000

yend=$( bc <<<"$ye - $timeint" )
for btime in $(seq $ys $timeint $yend)
do

     etime=$( bc <<<"$btime + $timeint" )
     couplingint=`echo "$timeint*1000.0" | bc -l`
     echo $btime $etime $timeint $couplingint 

     pismstarttime=$(date +%s)
     /bin/bash run_pism_step.sh ${SLURM_NTASKS} ${couplingint} ${it} >> pism/output/out_paleo_16km.out 
     pismendtime=$(date +%s) 

     srun --pack-group=1 python run_coupler.py -m "pism2vilma" -t ${timeint} -ys ${btime} -ye ${etime} -it ${it} --cpupertask=${CPUS_PER_TASK}

     vilmastarttime=$(date +%s)
     /bin/tcsh run_vilma3d_step.cli ${btime} ${etime} ${CPUS_PER_TASK}
     vilmaendtime=$(date +%s)

     srun --pack-group=1 python run_coupler.py -m "vilma2pism" -t ${timeint} -ys ${btime} -ye ${etime} -it ${it} --cpupertask=${CPUS_PER_TASK}
     endstep=$(date +%s)

     echo "STEP       $btime - $etime  kyr..."                         >> ./log-time.out
     echo "PISM       $(($pismendtime - $pismstarttime)) seconds..."   >> ./log-time.out
     echo "PISM2VILMA $(($vilmastarttime - $pismendtime))   seconds...">> ./log-time.out
     echo "VILMA      $(($vilmaendtime - $vilmastarttime))  seconds...">> ./log-time.out
     echo "VILMA2PISM $(($endstep - $vilmaendtime))   seconds..."      >> ./log-time.out
     echo "TOTAL STEP $(($endstep - $pismstarttime)) seconds..."       >> ./log-time.out
     echo                                                              >> ./log-time.out

     # zfill etime in output name #########
     if (( $(echo "$etime < 0.0" |bc -l) )); then
       e=$( bc <<<"$pad - $etime" )
       pismout="pism/results/pism2vilma-${e#?}-0.nc"
       pismres="pism/results/paleo-${e#?}-0.nc"
       vilmaout="out/vilma2pism-${e#?}-0.nc"
     else
       e=$( bc <<<"$pad*10 + $etime" )
       #echo $etime "pism2vilma${e#?}-0.nc"
       pismout="pism/results/pism2vilma${e#?}-0.nc"
       pismres="pism/results/paleo${e#?}-0.nc"
       vilmaout="out/vilma2pism${e#?}-0.nc"
     fi

     # check if outputfiles exist ##########
     if [ ! -f ${pismout} ]; then 
       echo "No PISM output file ${pismout}"
       break;
     fi
     if [ ! -f ${vilmaout} ]; then
       echo "No VILMA output file ${vilmaout}" 
       break;
     fi

done








