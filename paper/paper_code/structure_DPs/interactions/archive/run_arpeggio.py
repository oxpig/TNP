import os
import subprocess
from pathlib import Path

dataset_name = 'vhh_twist'
# path to NOT protonated structures - arpeggio will add hydrogens itself with flag 
structure_path = '/vols/bitbucket/gordon/tnp_data/vhh_twist_output/Raw_Model_Outputs/'
structures = [str(s).split('/')[-1] for s in list(Path(structure_path).glob('**/*.pdb'))][3:5]

# structure_path = '/data/localhost/gordon/TNP_Project/RESULTS/greiff_results/reduced_strucs/reduced_'+dataset_name+'/'
# structures = [str(s) for s in list(Path(structure_path).glob('*.pdb'))][0:1]

print(structures)
# structure_path = '/vols/bitbucket/gordon/tnp_data/vh_only_models/'+dataset_name+'/'
# structures = list(Path(structure_path).glob('*.pdb'))
# structures =  [s for s in structures if not s.name.endswith('_Annotated.pdb')][0:1] # remove ones annotated from TAP run

print(len(structures))

path_to_clean = '/data/localhost/gordon/TNP_Project/METRICS/greiff/gemma_scripts/arpeggio/clean_pdb.py'
path_to_arp = '/data/localhost/gordon/TNP_Project/METRICS/greiff/gemma_scripts/arpeggio/arpeggio.py'


def run_arpeggio(structure):
    #print(structure)
    # to fix from reruns
    bad = ['clean','hydrogenated']
    if any(b in str(structure) for b in bad):
        print('Ignoring:', str(structure))
    else:
        print('Calculating:', str(structure))
        # clean structure
        clean_cmd = 'python ' + path_to_clean + ' ' + structure_path + structure[:-4] + '/' + structure 
        #clean_cmd = 'python ' + path_to_clean + ' ' + str(structure)
        print(clean_cmd)
        subprocess.call(clean_cmd, shell = True) 
        # run arpeggio on clean structure
        run_cmd = 'python ' + path_to_arp + ' -v -wh -he ' + structure_path + structure[:-4] + '/' + structure[:-4] + '.clean.pdb' 
        #run_cmd = 'python ' + path_to_arp + ' -v -wh -he ' + str(structure_path) + str(structure).split('/')[-1][:-4] + '.clean.pdb'
        #run_cmd = 'python ' + path_to_arp + ' -v -he ' + str(structure)
        print(run_cmd)
        subprocess.call(run_cmd, shell = True) 


for structure in structures:
    run_arpeggio(structure)
