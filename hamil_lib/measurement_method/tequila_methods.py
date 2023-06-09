from typing import Tuple

from tequila.grouping.binary_rep import BinaryHamiltonian

from hamil_lib.measurement_method.measurement_method import MeasurementMethod


class TequilaMethods(MeasurementMethod):
    def __init__(self, method = "lf", options = None):
        super().__init__()
        self.method = method
        if options is None:
            options = {}
        self.options = options
        self.options["method"] = method

    def get_groups(self, H, options=None) -> Tuple[list, list]:
        Hbin = BinaryHamiltonian.init_from_qubit_hamiltonian(H)
        groups, ratios = Hbin.commuting_groups(options=self.options)
        if ratios[0] is not None:
            return groups, ratios
        for i in range(len(groups)):
            weight = 0.0
            for term in groups[i].binary_terms:
                weight += abs(term.coeff)
            ratios[i] = weight
        ratio_sum = sum(ratios)
        ratios = [ratio / ratio_sum for ratio in ratios]
        assert abs(sum(ratios) - 1.0) < 1e-6
        return groups, ratios