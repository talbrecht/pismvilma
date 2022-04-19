#!/bin/bash

export NN="$1"
export DT="$2"
export IT="$3"
IN=$((IT-1))

YS=-246000
YE=0

PISMR=pism/bin/pismr
SRUN="srun --pack-group=0 -n $NN"
OUTFOLDER=pism/results #-${IT}

mkdir -p $OUTFOLDER
mkdir -p pism/output

INPUTFOLDER=/p/tmp/albrecht/pism20/vilma/ressources/pism_input

ORIGFILE=$INPUTFOLDER/bedmap2_bheatflx_racmo_wessem_16km.nc
OKILLFILE=$INPUTFOLDER/okill_mask.nc
PICOFILE=$INPUTFOLDER/schmidtko14_edc-wdc_oceantemp_basins_response_fit_16km.nc
SLFILE=$INPUTFOLDER/imbrie06peltier15_sl.nc
DTFILE=$INPUTFOLDER/timeseries_edc-wdc_temp.nc

#INPUTFILE=$INPUTFOLDER/fit_16km_it0_tw.nc
#copy by run_coupler.py
INPUTFILE=pism/results/paleo_inp.nc


STRESS='-pik -sia_e 2.0 -sia_flow_law gpbld -limit_sia_diffusivity -stress_balance ssa+sia -ssa_method fd -ssa_flow_law gpbld -ssa_e 0.6'
SLIDE='-yield_stress mohr_coulomb -pseudo_plastic -pseudo_plastic_q 0.75 -pseudo_plastic_uthreshold 100.0 -hydrology null -tauc_slippery_grounding_lines'
CALV="-calving eigen_calving,thickness_calving -front_retreat_file $OKILLFILE -eigen_calving_K 1.0e17 -thickness_calving_threshold 75.0"
PICO="-ocean pico -ocean_pico_file $PICOFILE -gamma_T 1.0e-5 -overturning_coeff 0.8e6 -exclude_icerises -continental_shelf_depth -2500"
BEDDEF='-bed_def none'
SEALEV="-sea_level constant,delta_sl -ocean_delta_SL_file $SLFILE"
SMB="-atmosphere pik,delta_T,precip_scaling,elevation_change -atmosphere_pik era_interim_lon -atmosphere_pik_file $ORIGFILE -atmosphere_delta_T_file $DTFILE -atmosphere_precip_scaling_file $DTFILE -atmosphere_lapse_rate_file $ORIGFILE -temp_lapse_rate 8.2 -precip_adjustment scale -surface pdd"
TECH='-options_left -verbose 2 -backup_interval 3.0 -allow_extrapolation'
CONF="-config_override pism/conf/pism_config_override_pism1.2.nc"
TSERIES="-ts_file $OUTFOLDER/ts_paleo.nc -ts_times $YS:yearly:$YE"
RESULT="-o $OUTFOLDER/paleo.nc -o_size medium"
OUTPUT="$TSERIES $RESULT"

$SRUN $PISMR -i $INPUTFILE $STRESS $SLIDE $CALV $BEDDEF $SMB $TECH $CONF -y $DT $OUTPUT

