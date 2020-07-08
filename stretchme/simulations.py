from .stretching_tools import *


def simulate_traces(traces=1, p_prot=0.7, k_prot=200, p_dna=0, k_dna=None, position_blur=0.1, force_blur=1, l_dna=350,
                          l_prots=(25, 50, 100), rupture_forces=(10, 15), rupture_forces_blur=0.1,
                          force_range=(0.1, 20), relaxation=0.1):
    result = []
    for k in range(traces):
        result = simulate_single_trace(p_prot=p_prot, k_prot=k_prot, p_dna=p_dna, k_dna=k_dna, l_prots=l_prots,
                                       position_blur=position_blur, force_blur=force_blur, force_range=force_range,
                                       rupture_forces=rupture_forces, rupture_forces_blur=rupture_forces_blur,
                                       l_dna=l_dna, relaxation=relaxation)
    return result


def simulate_single_trace(p_prot=0.7, k_prot=200, p_dna=0, k_dna=0, position_blur=0.1, force_blur=1, l_dna=350,
                          l_prots=(25, 50, 100), rupture_forces=(10, 15), force_range=(0.1, 20), relaxation=0.1,
                          rupture_forces_blur=0.1):
    # TODO add validation of the input data
    # TODO make it more numpy-like
    # TODO we assume DNA relax in an instant. Is it right?
    result = pd.DataFrame(columns=['d', 'F'])
    real_rupture_forces = list(np.array(rupture_forces) + np.random.normal(0, rupture_forces_blur, len(l_prots)-1))
    real_rupture_forces.append(force_range[-1])
    d_prot_rupture = [l_prots[_] * (invert_wlc(real_rupture_forces[_], p_prot, k_prot) + relaxation)
                      for _ in range(len(l_prots) - 1)]
    start_forces = [force_range[0]] + [ewlc(d_prot_rupture[_], l_prots[_ + 1], p_prot, k_prot)
                                       for _ in range(len(d_prot_rupture))]
    part_forces = [[start_forces[_], real_rupture_forces[_]] for _ in range(len(real_rupture_forces))]

    for k in range(len(l_prots)):
        result = pd.concat([result, simulate_single_length(p_prot=p_prot, k_prot=k_prot, p_dna=p_dna, k_dna=k_dna,
                                                          l_prot=l_prots[k], l_dna=l_dna, position_blur=position_blur,
                                                          force_blur=force_blur, force_range=part_forces[k])])
    return result


def simulate_single_length(p_prot=0.7, k_prot=200, p_dna=0, k_dna=None, position_blur=0.1, force_blur=1, l_dna=350,
                          l_prot=25, force_range=(0.1, 20)):
    # TODO f_space cannot by < 0.1
    f_space = np.linspace(force_range[0], force_range[1], 10000)
    f_mixed = f_space + np.random.normal(0, force_blur, 10000)
    d_prot = l_prot * invert_wlc_np(f_space, p_prot, k_prot)
    if p_dna > 0:
        d_dna = l_dna * invert_wlc_np(f_space, p_dna, k_dna)
    else:
        d_dna = np.zeros(len(f_space))
    d_space = d_prot + d_dna
    d_mixed = d_space + np.random.normal(0, position_blur, 10000)
    trace = pd.DataFrame({'d': d_mixed, 'F': f_mixed})
    return trace
