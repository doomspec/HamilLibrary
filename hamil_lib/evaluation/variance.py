from openfermion import bravyi_kitaev
from tequila import QubitHamiltonian

from hamil_lib.evaluation.evaluate import get_variance_on_average_state
from hamil_lib.measurement_method.l1_method import L1Method
from hamil_lib.measurement_method.measurement_method import MeasurementMethod
from hamil_lib.measurement_method.tequila_methods import TequilaMethods


def evaluate_avg_var_on_small_mols(method: MeasurementMethod, transformation, mol_iterator):
    vars = {}
    for op, record in mol_iterator():
        H = QubitHamiltonian.from_openfermion(transformation(op))
        var = get_variance_on_average_state(H, method)
        vars[record["name"]] = var
        print(record["name"]+":"+str(var))
    return vars

if __name__ == '__main__1':
    from hamil_lib.molecule.small_mols import iter_small_mols

    evaluate_avg_var_on_small_mols(L1Method(), bravyi_kitaev,
                                   iter_small_mols)
    evaluate_avg_var_on_small_mols(TequilaMethods("si", {"condition": "qwc"}), bravyi_kitaev,
                                   iter_small_mols)
    evaluate_avg_var_on_small_mols(TequilaMethods("si", {"condition": "fc"}), bravyi_kitaev, iter_small_mols)

"""
H2_6-31g:130.72265215675762
LiH_sto-3g:79.03375274163545
LiH_6-31g:1492.9423865157069
H2O_6-31g:12779.819222079866
H2_6-31g:22.531463303985205
LiH_sto-3g:9.343076656897285
LiH_6-31g:164.188711177738
H2O_6-31g:2808.5659241747394
H2_6-31g:8.548648733032492
LiH_sto-3g:4.614736866420189
LiH_6-31g:39.90999149264847
H2O_6-31g:230.1369499426702
"""

if __name__ == '__main__':
    from hamil_lib.molecule.small_mols import iter_small_mols
    evaluate_avg_var_on_small_mols(TequilaMethods("ics", {"condition": "fc"}), bravyi_kitaev, iter_small_mols)