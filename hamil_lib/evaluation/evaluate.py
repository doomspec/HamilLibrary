from tequila import QubitHamiltonian
from tequila.grouping.binary_rep import BinaryPauliString, BinaryHamiltonian

from hamil_lib.measurement_method.measurement_method import MeasurementMethod
from hamil_lib.evaluation.wfn_utils import get_cov_dict


def get_cover_probs(groups, ratios):
    joint_prob_dict = {}
    prob_dict = {}
    for idx_g, part in enumerate(groups):
        op: QubitHamiltonian
        opbin = part
        part_prob = ratios[idx_g]
        group_terms = opbin.binary_terms
        for idx1, term1 in enumerate(group_terms):
            pw1 = BinaryPauliString(term1.get_binary(), 1.0).binary_tuple()
            if pw1 in prob_dict.keys():
                prob_dict[pw1] += part_prob
            else:
                prob_dict[pw1] = part_prob
            for idx2, term2 in enumerate(group_terms[idx1 + 1:]):
                pw2 = BinaryPauliString(term2.get_binary(), 1.0).binary_tuple()
                # assert term1.commute(term2)
                pw_pair = (pw1, pw2)
                r_pw_pair = (pw2, pw1)
                if pw_pair not in joint_prob_dict.keys():
                    joint_prob_dict[pw_pair] = part_prob
                    # assert r_pw_pair not in joint_prob_dict.keys()
                    joint_prob_dict[r_pw_pair] = part_prob
                else:
                    joint_prob_dict[pw_pair] += part_prob
                    joint_prob_dict[r_pw_pair] += part_prob
    return joint_prob_dict, prob_dict

def get_coeff_dict(Hbin):
    coeff_dict = {}
    Hbin_terms = Hbin.binary_terms[1:]
    for idx, term in enumerate(Hbin_terms):
        pw = BinaryPauliString(term.get_binary(), 1.0).binary_tuple()
        coeff_dict[pw] = Hbin_terms[idx].coeff
    return coeff_dict

def calc_variance(Hbin, joint_prob_dict, prob_dict, cov_dict):
    var_coeff = 0.0
    Hbin_terms = Hbin.binary_terms[1:]
    seen_pairs = set()
    for idx1, term1 in enumerate(Hbin_terms):
        for idx2, term2 in enumerate(Hbin_terms[idx1+1:]):
            pw1 = BinaryPauliString(term1.get_binary(), 1.0).binary_tuple()
            pw2 = BinaryPauliString(term2.get_binary(), 1.0).binary_tuple()
            pw_pair = (pw1, pw2)
            seen_pairs.add(pw_pair)
            assert (pw2, pw1) not in seen_pairs
            assert pw1 != pw2
            joint_prob = joint_prob_dict.get(pw_pair, 0.0)
            coeff_prod = term1.coeff * term2.coeff
            cov = cov_dict.get(pw_pair, 0.0)
            if cov == 0.0 or joint_prob == 0.0:
                continue
            prob_prod = prob_dict.get(pw1, 0.0) * prob_dict.get(pw2, 0.0)
            assert prob_prod != 0.0
            prob = joint_prob / prob_prod
            inc = coeff_prod * cov * prob
            # Multiply by 2 because we only considered half of the covariance matrix
            var_coeff += 2.0 * inc

    for idx, term in enumerate(Hbin_terms):
        pw = BinaryPauliString(term.get_binary(), 1.0).binary_tuple()
        coeff = term.coeff
        pw_pair = (pw, pw)
        cov = cov_dict[pw_pair]
        prob = prob_dict.get(pw, 0.0)
        coeff2 = coeff ** 2
        var_coeff += coeff2 * cov / prob

    return var_coeff

def get_variance_on_ground_state(H: QubitHamiltonian, measurement_method: MeasurementMethod):
    Hbin = BinaryHamiltonian.init_from_qubit_hamiltonian(H)
    assert max(Hbin.binary_terms[0].binary_tuple()) == 0.0
    groups, ratios = measurement_method.get_groups(H)
    joint_prob_dict, prob_dict = get_cover_probs(groups, ratios)
    cov_dict = get_cov_dict(Hbin)
    variance = calc_variance(Hbin, joint_prob_dict, prob_dict, cov_dict)
    return variance

def calc_average_variance(Hbin,prob_dict):

    Hbin_terms = Hbin.binary_terms[1:]
    diagonal_var = 2**Hbin.n_qubit
    diagonal_var = diagonal_var / (diagonal_var + 1)
    var_coeff = 0.0
    for idx, term in enumerate(Hbin_terms):
        pw = BinaryPauliString(term.get_binary(), 1.0).binary_tuple()
        coeff = term.coeff
        prob = prob_dict.get(pw, 0.0)
        coeff2 = coeff ** 2
        var_coeff += coeff2 * diagonal_var / prob

    return var_coeff

def get_variance_on_average_state(H: QubitHamiltonian, measurement_method: MeasurementMethod):
    Hbin = BinaryHamiltonian.init_from_qubit_hamiltonian(H)
    assert max(Hbin.binary_terms[0].binary_tuple()) == 0.0
    groups, ratios = measurement_method.get_groups(H)
    joint_prob_dict, prob_dict = get_cover_probs(groups, ratios)
    variance = calc_average_variance(Hbin, prob_dict)
    return variance

