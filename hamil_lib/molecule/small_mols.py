import os

from hamil_lib_helper.hamil import read_hamil_record


def iter_small_mols():
    """Iterate over all small molecules in the database."""
    _small_mols = [
        "H2_6-31g",
        "LiH_sto-3g",
    ]
    mol_folder_path = os.path.dirname(os.path.abspath(__file__))+"/output"
    for name in _small_mols:
        op, record = read_hamil_record(mol_folder_path + '/' + name + '.hamil')
        yield op, record


if __name__ == '__main__':
    for op, record in iter_small_mols():
        print(record)