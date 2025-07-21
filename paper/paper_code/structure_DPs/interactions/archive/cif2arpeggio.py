from multiprocessing import Pool
import subprocess
import os
import time

cpus = 80

structure_path = '/storage/evagsm/nobackup/crystal_dataset/structures/AF2_multi_L_cif/'
output_path = '/storage/evagsm/nobackup/crystal_dataset/structures/AF2_multi_L_arpeggio/'
#l1 = os.listdir(structure_path) 
#l2 = os.listdir(output_path) 
#structures = [x for x in l1 if x.split('.')[0] + '.json' not in l2]

structures = os.listdir(structure_path)

#print(len(l1))
#print(len(l2))
print(len(structures))


cmds = ['OUTPUT_PATH=/storage/evagsm/abodybuilder/arpeggio/', 'cd $OUTPUT_PATH']
subprocess.call(";".join(cmds), shell = True)    
    
def run(structure):
    print(structure)
    cmd = ('arpeggio {}'.format(structure_path + structure))
    subprocess.call(cmd, shell = True) 
      
if __name__ == '__main__':
    pool = Pool(processes = cpus)
    results = pool.map(run, structures)
    pool.close()
    pool.join()
