from copy import deepcopy
from typing import Tuple

from openfermion import QubitOperator
from tequila.grouping.binary_rep import BinaryHamiltonian

from hamil_lib.measurement_method.measurement_method import MeasurementMethod


class L1Method(MeasurementMethod):
    def __init__(self):
        super().__init__()

    def get_groups(self, H: QubitOperator, options=None) -> Tuple[list, list]:
        assert () not in H.terms
        groups = []
        ratios = []
        for ps, coeff in H.terms.items():
            groups.append([ps])
            ratios.append(abs(coeff))
        ratio_sum = sum(ratios)
        ratios = [ratio / ratio_sum for ratio in ratios]
        return groups, ratios