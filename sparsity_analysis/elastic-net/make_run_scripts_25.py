#!/usr/bin/env python

'''make a run script for each subset and output a qsub file'''

qsubfile = open('../qsub.txt','w')
prescript = '25_GTEx_OTD_TS_CV_elasticNet'
tissues = [line.rstrip('\n') for line in open('/group/im-lab/nas40t2/hwheeler/cross-tissue/ten.tissue.list')]

#NOTE: walltime of 24hr sufficient for glmnet

for i in range(1,23):
    newi = str(i)
    for alpha in [0.05, 0.5, 0.95, 1]: 
        outfilename = 'run_' + prescript + '_' + newi + '.' + str(alpha) + '.sh'
        outfile = open(outfilename,'w')
        output = '''#!/bin/bash
#PBS -N R.glmnet.''' + newi + '_' + str(alpha) + '''\n#PBS -S /bin/bash
#PBS -l walltime=240:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=12gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

module load R\n\n'''

        outfile.write(output)
        for tis in tissues[1:]:
            command = 'time R --no-save < ''' + prescript + '.r --args ' + newi + ' ' + str(alpha) + ' \"' + tis + '\"\n'
            outfile.write(command)
        
        qsubfile.write('qsub run_scripts/' + outfilename + '\nsleep 3\n')


