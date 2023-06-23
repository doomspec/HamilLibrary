from openfermion import bravyi_kitaev
from tequila import QubitHamiltonian

from hamil_lib.evaluation.evaluate import get_variance_on_average_state
from hamil_lib.measurement_method.l1_method import L1Method
from hamil_lib.measurement_method.measurement_method import MeasurementMethod
from hamil_lib.measurement_method.tequila_methods import TequilaMethods


def evaluate_on_small_mols(method: MeasurementMethod, transformation, mol_iterator):
    vars = {}
    for op, record in mol_iterator():
        H = QubitHamiltonian.from_openfermion(transformation(op))
        var = get_variance_on_average_state(H, method)
        vars[record["name"]] = var
        print(record["name"]+":"+str(var))
    return vars

if __name__ == '__main__':
    from hamil_lib.molecule.small_mols import iter_small_mols

    evaluate_on_small_mols(L1Method(), bravyi_kitaev,
                           iter_small_mols)
    evaluate_on_small_mols(TequilaMethods("si", {"condition": "qwc"}), bravyi_kitaev,
                           iter_small_mols)
    evaluate_on_small_mols(TequilaMethods("si", {"condition": "fc"}), bravyi_kitaev, iter_small_mols)
