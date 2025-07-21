import os
import subprocess
from pathlib import Path
from multiprocessing import Pool

cpus = 23

# running on hydrogenated structures
# from greiff paper:
# Finally, we used Arpeggio (version 1.4.1) to calculate interatomic interactions 
#   after converting antibody hydrogenated structures to cif format
dataset_name = 'vhh_tsd'
structure_path = '/vols/bitbucket/gordon/tnp_data/reduced_strucs/reduced_'+dataset_name
structures = [s for s in os.listdir(structure_path)]
output_path = '/vols/bitbucket/gordon/tnp_data/arpeggio/'

print('TOTAL STRUCTURES:', len(structures))

def run_arpeggio(structure):

    # to fix from reruns
    try:
        print('Calculating:', str(structure))
        
        # convert from pdb to mmcif format which newer version of arpeggio depends on 
        convert_cmd = 'gemmi convert ' + structure_path + '/' + structure + ' ' +  structure_path + '/' + structure[:-4] + '.cif'
        #print(convert_cmd)
        subprocess.call(convert_cmd, shell = True) 
        
        # run arpeggio on cif file
        run_cmd = 'pdbe-arpeggio -m ' + structure_path + '/' + structure[:-4] + '.cif -o ' + output_path + dataset_name + '/'
        #print(run_cmd)
        subprocess.call(run_cmd, shell = True) 
    
    except:
        print('FAILED:', structure)

# for structure in structures:
#     run_arpeggio(structure)

if __name__ == '__main__':
    pool = Pool(processes = cpus)
    results = pool.map(run_arpeggio, structures)
    pool.close()
    pool.join()