from .stretching_tools import *


def simulate_traces(traces=1, p_prot=0.7, k_prot=200, p_dna=0, k_dna=None, position_blur=0.1, force_blur=0.1,
                    l_prots=(50, 100), rupture_forces=(10), rupture_forces_blur=(0.1), force_range=(0, 20)):
    result = []
    for k in range(traces):
        result = simulate_single_trace(p_prot=p_prot, k_prot=k_prot, p_dna=p_dna, k_dna=k_dna, l_prots=l_prots,
                                       position_blur=position_blur, force_blur=force_blur, force_range=force_range,
                                       rupture_forces=rupture_forces, rupture_forces_blur=rupture_forces_blur)
    return result


def simulate_single_trace(p_prot=0.7, k_prot=200, p_dna=0, k_dna=None, position_blur=0.1, force_blur=0.1,
                          l_prots=(50, 100), rupture_forces=(10), rupture_forces_blur=(0.1), force_range=(0, 20)):
    f_space = np.linspace(force_range[0], force_range[1], 1000)
    f_mixed = f_space + np.random.normal(0, force_blur, 1000)
    d_prot = invert_wlc_np(f_space, p_prot, k_prot)
    if p_dna > 0:
        d_dna = invert_wlc_np(f_space, p_dna, k_dna)
    else:
        d_dna = np.zeros(len(f_space))
    d_space = d_prot + d_dna
    d_mixed = d_space + np.random.normal(0, position_blur, 1000)
    trace = pd.DataFrame({'d': d_mixed, 'F': f_mixed})
    return trace
