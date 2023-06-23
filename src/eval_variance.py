
from tequila.grouping.binary_rep import BinaryHamiltonian
import tequila as tq
from measurement_method.evaluate import get_variance_on_ground_state
from measurement_method.l1_method import L1Method
from measurement_method.tequila_methods import TequilaMethods

basis_set = 'sto-3g'
molecule = tq.quantumchemistry.Molecule(geometry="./molecule/geometry/LiH.xyz", basis_set=basis_set, transformation="bravyi-kitaev")

H = molecule.make_hamiltonian()
n_qubits = H.n_qubits
Hbin = BinaryHamiltonian.init_from_qubit_hamiltonian(H)

var = get_variance_on_ground_state(H, TequilaMethods("si", {"condition": "qwc"}))
var2 = get_variance_on_ground_state(H, L1Method())
var3 = get_variance_on_ground_state(H, TequilaMethods("si"))

weight = 0.0
for term in Hbin.binary_terms[1:]:
    weight += abs(term.coeff)

print(var)
print(var2)
print(var3)
print(weight**2)