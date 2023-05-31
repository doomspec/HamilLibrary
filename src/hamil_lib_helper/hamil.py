import os
import pickle
import openfermion
from openfermion import SymbolicOperator, QubitOperator, FermionOperator


def make_hamil_record(op: SymbolicOperator, type, n_qubit, basis):
    assert type in ["fermion", "qubit"]
    assert n_qubit > 0
    return {
        "terms": op.terms,
        "type": type,
        "n_qubit": n_qubit,
        "basis": basis
    }


def read_hamil_record(file_path):
    with open(file_path, "rb") as f:
        record = pickle.load(f)
    op = None
    if record["type"] == "fermion":
        op = FermionOperator()
    elif record["type"] == "qubit":
        op = QubitOperator()
    op.terms = record["terms"]
    return op, record


def save_fermion_hamil_for_mol(mol, basis_set, emitting_file):
    fop = make_fermionic_hamiltonian(mol)
    record = make_hamil_record(fop, "fermion", mol.n_orbitals * 2, basis_set)
    file_path = os.path.abspath(emitting_file)
    file_name = ".".join(os.path.basename(file_path).split(".")[:-1])
    dir_path = os.path.dirname(file_path)
    with open(dir_path + '/output/' + file_name + '.hamil', 'wb') as f:
        pickle.dump(record, f)


def make_fermionic_hamiltonian(mol) -> SymbolicOperator:
    of_molecule = mol.make_molecule()
    fop = of_molecule.get_molecular_hamiltonian()
    fop = openfermion.transforms.get_fermion_operator(fop)
    return fop
