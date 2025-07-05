from openfermion import bravyi_kitaev

from hamil_lib.evaluation.variance import evaluate_avg_var_on_small_mols
from hamil_lib.measurement_method.l1_method import L1Method
from hamil_lib.measurement_method.tequila_methods import TequilaMethods

if __name__ == '__main__':
    from hamil_lib.molecule.small_mols import iter_small_mols

    #evaluate_avg_var_on_small_mols(L1Method(), bravyi_kitaev,
    #                               iter_small_mols)
    evaluate_avg_var_on_small_mols(TequilaMethods("si", {"condition": "qwc"}), bravyi_kitaev,
                                   iter_small_mols)
    evaluate_avg_var_on_small_mols(TequilaMethods("si", {"condition": "fc"}), bravyi_kitaev, iter_small_mols)
