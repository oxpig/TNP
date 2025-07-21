import os
from Bio.PDB import PDBParser, MMCIFIO

from multiprocessing import Pool


cpus = 80


structure_path = '/storage/evagsm/nobackup/crystal_dataset/structures/AF2_multi_L_h/'
output_path = '/storage/evagsm/nobackup/crystal_dataset/structures/AF2_multi_L_cif/'

structures = os.listdir(structure_path)

def pdb2cif(pdb, structure_path = structure_path):
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('', structure_path + pdb) 
        io = MMCIFIO()
        io.set_structure(structure)
        io.save("{}.cif".format(output_path + pdb))
        print('{} done'.format(pdb))
    except Exception:
        print('{} failed'.format(pdb))


if __name__ == '__main__':
    pool = Pool(processes = cpus)
    results = pool.map(pdb2cif, structures)
    pool.close()
    pool.join()
