import os
from os.path import exists,isdir
from molSimplify.job_manager.tools import get_total_queue_usage, list_active_jobs, call_bash


def read_Orca(filepath):

    #Check if optimization has converged normally (not inconvergence or out of time)
    #If so,return the final energy: energy + external potential

    if exists(filepath) == False:
        return 'Failed_nofile'

    lines = open(filepath,'r').readlines()
    cond_1 = False #'Optimization run done' is found in file
    cond_2 = False #'Orca terminated normally' is found (as second-last line)
    cond_3 = True # 'The optimization did not converge' is found in file

    Esp,Eep = [],[]
    for line in lines:

        if 'OPTIMIZATION RUN DONE' in line:
            cond_1 = True

        elif 'ORCA TERMINATED NORMALLY' in line:
            cond_2 = True

        elif 'The optimization did not converge' in line:
            cond_3 = False

        elif 'FINAL SINGLE POINT ENERGY' in line:
            Esp.append(float(line.split()[-1]))

        elif 'External Potential' in line:
            Eep.append(float(line.split()[3]))

    if cond_1 and cond_2 and cond_3:
        if len(Eep) != 0:#might be empty for 0nN
            return Esp[-1] + Eep[-1]
        else:
            return Esp[-1]
    else:
        return 'Failed'

#Obtain number of jobs to be submitted this call
max_subs = 200
remaining_subs = max_subs - get_total_queue_usage()

all_jobs = []
for i in os.listdir():
    if isdir(i):
        all_jobs.append(i)
active_jobs = list_active_jobs()
inactive_jobs = [x for x in all_jobs if x not in active_jobs]

tosubmit_jobs = []
for job in inactive_jobs:
    outpath = './' + job + '/' + job + '.out'
    if exists(outpath) == False:
        tosubmit_jobs.append(job)
    elif type(read_Orca(outpath)) == str:
        tosubmit_jobs.append(job)

submission = remaining_subs
idx = 0
while submission > 0:
    job = tosubmit_jobs[idx]
    os.chdir('./'+job)
    call_bash('sbatch ' + job + '_jobscript',error = True)
    os.chdir('..')
    idx += 1
    submission -= 1

print(str(remaining_subs) + ' jobs submitted')
