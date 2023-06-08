from copy import deepcopy
from typing import Tuple

from tequila.grouping.binary_rep import BinaryHamiltonian, BinaryPauliString

from measurement_method.measurement_method import MeasurementMethod


class L1Method(MeasurementMethod):
    def __init__(self):
        super().__init__()

    def get_groups(self, H, options=None) -> Tuple[list, list]:
        Hbin = BinaryHamiltonian.init_from_qubit_hamiltonian(H)
        groups = []
        ratios = []
        for term in Hbin.binary_terms[1:]:
            groups.append(BinaryHamiltonian(binary_terms=[deepcopy(term)]))
            ratios.append(abs(term.coeff))
        ratio_sum = sum(ratios)
        ratios = [ratio / ratio_sum for ratio in ratios]
        assert abs(sum(ratios) - 1.0) < 1e-6
        return groups, ratios