#from multiprocessing import Pool
import subprocess
import os
#import time

import os
import subprocess
from pathlib import Path

dataset_name = 'vhh_twist'

# structure_path = '/storage/evagsm/nobackup/crystal_dataset/structures/AF2_multi/L/'
# output_path = '/storage/evagsm/nobackup/crystal_dataset/structures/AF2_multi_L_h/'
# structures = os.listdir(structure_path)

structure_path = '/vols/bitbucket/gordon/tnp_data/vhh_twist_output/Raw_Model_Outputs/'
output_path = '/data/localhost/gordon/TNP_Project/RESULTS/greiff_results/reduced_strucs/reduced_'+dataset_name+'/'
structures = list(Path(structure_path).glob('**/*.pdb'))
#structures =  [s for s in structures if not s.name.endswith('_Annotated.pdb')] # remove ones annotated from TAP run

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




