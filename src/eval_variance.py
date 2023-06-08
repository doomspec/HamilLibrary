
from tequila.grouping.binary_rep import BinaryHamiltonian
import tequila as tq
from measurement_method.evaluate import get_variance_on_ground_state
from measurement_method.tequila_methods import TequilaMethods

basis_set = 'sto-3g'
molecule = tq.quantumchemistry.Molecule(geometry="./molecule/geometry/LiH_short.xyz", basis_set=basis_set, transformation="bravyi-kitaev")

H = molecule.make_hamiltonian()
n_qubits = H.n_qubits
Hbin = BinaryHamiltonian.init_from_qubit_hamiltonian(H)

var = get_variance_on_ground_state(H, TequilaMethods("lf"))
weight = 0.0
for term in Hbin.binary_terms:
    weight += abs(term.coeff)
print(var)
print(weight**2)