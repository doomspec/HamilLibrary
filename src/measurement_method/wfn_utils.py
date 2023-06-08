import os
import pickle
import hashlib
from openfermion import qubit_operator_sparse, get_ground_state, get_sparse_operator, count_qubits
from tequila import QubitHamiltonian
from tequila.grouping.binary_rep import BinaryPauliString, BinaryHamiltonian
import tequila as tq


def get_cov_dict_from_psi(Hbin, psi):
    '''
    Get <P_1P_2>-<P_1><P_2> for all P_1, P_2 in H
    '''
    terms = Hbin.binary_terms
    cov_dict = {}
    wfn0 = tq.QubitWaveFunction(psi)
    for idx, term1 in enumerate(terms):
        for term2 in terms[idx:]:
            pw1 = BinaryPauliString(term1.get_binary(), 1.0)
            pw2 = BinaryPauliString(term2.get_binary(), 1.0)
            op1 = QubitHamiltonian.from_paulistrings(pw1.to_pauli_strings())
            op2 = QubitHamiltonian.from_paulistrings(pw2.to_pauli_strings())
            if pw1.commute(pw2):
                prod_op = op1 * op2
                cov_dict[(term1.binary_tuple(), term2.binary_tuple())] = wfn0.inner(prod_op(wfn0)) - wfn0.inner(op1(wfn0)) * wfn0.inner(op2(wfn0))
    return cov_dict


def get_cov_dict(Hbin: BinaryHamiltonian):
    H = Hbin.to_qubit_hamiltonian()
    H_hash = hashlib.md5(str(H).encode("utf-8")).hexdigest()
    cache_path = os.path.dirname(os.path.abspath(__file__))+"/cache_cov_dict_"+H_hash+".pkl"
    if os.path.isfile(cache_path):
        with open(cache_path, 'rb') as file:
            cov_dict = pickle.load(file)
        return cov_dict
    sparse_op = qubit_operator_sparse(Hbin.to_qubit_hamiltonian().qubit_operator)
    eigv, psi = get_ground_state(sparse_op)

    cov_dict = get_cov_dict_from_psi(Hbin, psi)

    with open(cache_path, 'wb') as file:
        pickle.dump(cov_dict, file)

    return cov_dict