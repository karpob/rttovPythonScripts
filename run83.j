#!/bin/csh -fx
# ------------------------------
#SBATCH -A g0613
#SBATCH --export=NONE
#
#PBS -N run83
#PBS -o run83.log.o%j
#SBATCH --ntasks=1
##SBATCH --ntasks-per-node=16
#SBATCH --constraint=hasw
##SBATCH --qos=debug
##PBS -l walltime=1:00:00
##SBATCH --partition=preops
##SBATCH --qos=dastest
##SBATCH --qos=obsdev
#SBATCH --qos=debug
#PBS -l walltime=0:50:00
##PBS -l mem=4gb
#PBS -S /bin/csh
#PBS -j eo
#BSUB -J m2m1c
#BSUB -n 384
#BSUB -W 5:00
#BSUB -o m2m1c.log.o%J
#BSUB -e m2m1c.log.o%J
# ------------------------------
source  /discover/nobackup/bkarpowi/github/rttovPythonScripts/sourceMe 
which python3 
python3 /discover/nobackup/bkarpowi/github/rttovPythonScripts/iasiCrisAirs_inBufrSubset_allStandard.py
