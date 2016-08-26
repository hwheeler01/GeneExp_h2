#!/bin/bash
#PBS -N ldak
#PBS -S /bin/bash
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=4gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

##THIS IS A TEST SCRIPT
##SEE run_scripts/make_run_scripts_01.py & run_scripts/run_01... FOR USED SCRIPTS

for i in {1..2}
do
    ./ldak.4.9 --calc-weights examples/sections --bfile examples/test --section $i
done

./ldak.4.9 --join-weights examples/sections
