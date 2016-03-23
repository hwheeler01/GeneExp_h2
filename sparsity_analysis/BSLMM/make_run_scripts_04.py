#!/usr/bin/env python

'''make a run script for each subset and output a qsub file'''
import string
qsubfile = open('../qsub.txt','w')
prescript = '04_gtex_tissue-wide_bslmm_A-P'
tissues = [line.rstrip('\n') for line in open('/group/im-lab/nas40t2/hwheeler/cross-tissue/nine.tissue.list')]

for i in range(1,23): 
    newi = str(i)
    for alpha in list(string.ascii_uppercase)[0:16]: #divide each chr into 16 parts
        outfilename = 'run_' + prescript + '_' + newi + '.' + str(alpha) + '.sh'
        outfile = open(outfilename,'w')
        output = '''#!/bin/bash
#PBS -N R.bslmm.''' + newi + '.' + str(alpha) + '''\n#PBS -S /bin/bash
#PBS -l walltime=480:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

module load R/3.2.0
module load gemma/0.94\n\n'''

        outfile.write(output)
        for tis in tissues:
            command = 'time R --no-save < ''' + prescript + '.r --args ' + newi + ' ' + str(alpha) + ' \"' + tis + '\"\n'
            outfile.write(command)

        qsubfile.write('qsub run_scripts/' + outfilename + '\nsleep 3\n')


