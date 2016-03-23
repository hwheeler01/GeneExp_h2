#!/usr/bin/env python

'''make a run script for each subset and output a qsub file'''
import string
qsubfile = open('../qsub.txt','w')
prescript = '01_dgn_bslmm'

for i in range(1,23):
    newi = str(i)
    for alpha in list(string.ascii_uppercase)[0:8]: #divide each chr into 8 parts
        outfilename = 'run_' + prescript + '_' + newi + '.' + str(alpha) + '.sh'
        outfile = open(outfilename,'w')
        output = '''#!/bin/bash
#PBS -N R.bslmm.''' + newi + '.' + str(alpha) + '''\n#PBS -S /bin/bash
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

module load R
module load gemma/0.94

time R --no-save < ''' + prescript + '.r --args ' + newi + ' ' + str(alpha) + '\n'

        outfile.write(output)
        qsubfile.write('qsub run_scripts/' + outfilename + '\nsleep 3\n')


