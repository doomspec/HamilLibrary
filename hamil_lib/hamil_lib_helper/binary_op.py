from openfermion import QubitOperator, count_qubits


class BinaryHamiltonian:
    def __init__(self, terms: dict, n_qubits: int, ignore_const=False):
        '''
        Initiate from a list of Binary Pauli Strings
        '''

        self.n_qubit = n_qubits



def terms_to_binary_terms(terms: dict, n_qubit):
    bin_terms = []
    coeff = []
    for op_string, coeff in terms.items():
        coeff.append(coeff)
        bin_term = [0] * (2 * n_qubit)
        bin_terms.append(bin_term)
        for idx, pauli in op_string:
            if pauli == "X":
                bin_term[idx] = 1
            elif pauli == "Y":
                bin_term[idx] = 1
                bin_term[n_qubit + idx] = 1
            elif pauli == "Z":
                bin_term[n_qubit + idx] = 1


if __name__ == "__main__":
    op = QubitOperator(((0, "X"), (1, "Y"), (2, "Z")))
    op = QubitOperator("X0 Y1 Z2") + QubitOperator("X0 Y2 Z4")
    n_qubit = count_qubits(op)
    print(op.terms)
