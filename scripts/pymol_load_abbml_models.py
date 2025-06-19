"""Execution $ pymol pymol_load_abbml_models.py -- folderpath_with_pdbs"""
import os, sys
from pymol import cmd

workdir = sys.argv[1]

pdb_id_list = [item for item in os.listdir(workdir) if os.path.isdir(os.path.join(workdir,item))]
pdb_path_suffix = "my_query/output/my_query_rank1_imgt_scheme.pdb"
pdb_path_list = [os.path.join(*[workdir, pdb_id, pdb_path_suffix]) for pdb_id in pdb_id_list]

for i in range(len(pdb_id_list)):
    pdb_path = pdb_path_list[i]
    pdb_id = pdb_id_list[i]
    try:
        cmd.load(pdb_path, pdb_id)
    except:
        print("Model failed to load ", pdb_id)
