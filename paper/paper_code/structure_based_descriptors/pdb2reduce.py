#from multiprocessing import Pool
import subprocess
import os
#import time

import os
import subprocess
from pathlib import Path

dataset_name = ''

structure_path = '/models/'
output_path = '/data/localhost/gordon/TNP_Project/RESULTS/greiff_results/reduced_strucs/reduced_'+dataset_name+'/'
structures = list(Path(structure_path).glob('**/*.pdb')) 

print(len(structures))


def run_reduce(structure):

    # print(type(structure))
    # c
    # else:
    print(structure)
    #cmd = 'reduce -FLIP {} > {}'.format(structure_path + structure[:-4] + '/' + structure, output_path + structure[:-4])
    cmd = 'reduce -FLIP {} > {}'.format(structure, output_path + str(structure).split('/')[-1][:-4])
    #print(cmd)
    subprocess.call(cmd, shell = True) 

for structure in structures:
    run_reduce(structure)




