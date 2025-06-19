import re
from Bio.PDB import PDBParser, PDBIO, Select

_hydrogen = re.compile("[123 ]*H.*")

class RemoveHydrogen(Select):
    def accept_atom(self, atom):
        """Verify if atoms are not Hydrogen."""
        # atoms - get rid of hydrogens
        name = atom.get_id()
        return 0 if _hydrogen.match(name) else 1

def pdb_remove_hydrogens(input_pdb):
    try:
        pdb = PDBParser().get_structure("antibody", input_pdb)
        io = PDBIO()
        io.set_structure(pdb)

        output_model = f'{input_pdb[:-4]}.pdb'
        io.save(output_model, RemoveHydrogen())

    except Exception as e:
        print(e)

if __name__ == "__main__":
    import sys
    input_pdb = sys.argv[1]
    pdb_remove_hydrogens(input_pdb)
