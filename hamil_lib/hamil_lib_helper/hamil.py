import os
import pickle
import openfermion
from openfermion import SymbolicOperator, QubitOperator, FermionOperator
from openfermion import count_qubits


def make_hamil_record(op: SymbolicOperator, type, n_qubit, **kwargs):
    assert type in ["fermion", "qubit"]
    assert n_qubit > 0
    res = {
        "terms": op.terms,
        "type": type,
        "n_qubit": n_qubit,
    }
    for key, value in kwargs.items():
        assert key in ["basis"]
        res[key] = value
    return res

def read_hamil_record(file_path) -> [SymbolicOperator, dict]:
    with open(file_path, "rb") as f:
        record = pickle.load(f)
    op = None
    if record["type"] == "fermion":
        op = FermionOperator()
    elif record["type"] == "qubit":
        op = QubitOperator()
    op.terms = record["terms"]
    record["name"] = os.path.basename(file_path).split(".")[0]
    del record["terms"]
    return op, record


def save_fermion_hamil(mol, basis_set, emitting_file):
    fop = make_fermionic_hamiltonian(mol)
    record = make_hamil_record(fop, "fermion", mol.n_orbitals * 2, basis=basis_set)
    file_path = os.path.abspath(emitting_file)
    file_name = ".".join(os.path.basename(file_path).split(".")[:-1])
    dir_path = os.path.dirname(file_path)
    with open(dir_path + '/output/' + file_name + '.hamil', 'wb') as f:
        pickle.dump(record, f)


def save_qubit_hamil(op: QubitOperator, emitting_file):
    record = make_hamil_record(op, "qubit", count_qubits(op))
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
