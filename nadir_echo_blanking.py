"""
# author Fabian Glatz
##################################################################
# Description:
# This two functions shows a first idea how the nadir echo can be detected and removed during dual-focus post-processing.
# These functions are performed after matched filtering with a filter adapted to the nadir echo (nadir echo focused)
##################################################################
"""
import matplotlib.pyplot as plt
import numpy as np



def nadir_echo_blanking(u_foc_nadir,u1_foc_nadir, u2_foc_nadir, n_average, blanking_dB_rel_max_amp):
    blanking_value_rel_max=10**(blanking_dB_rel_max_amp/20)
    size=u2_foc_nadir.shape
    n_times=np.int(np.round(size[1]/n_average))
    num_blanked_lines = np.zeros(n_times)
    for k in range(0, n_times):
        print(k)
        part=u2_foc_nadir[:, k * n_average:(k + 1) * n_average]
        size_part=part.shape
        nadir_shape = np.abs(part).sum(axis=1)/size_part[1]
        mask = np.where(nadir_shape > nadir_shape.max() / blanking_value_rel_max, 0, 1)
        blanked_lines=np.array(np.where(mask == 0))
        start=blanked_lines[0,0]
        stop=blanked_lines[0,-1]+1
        num_blanked_lines[k]= np.count_nonzero(mask == 0)
        u_foc_nadir[start:stop,k * n_average:(k + 1) * n_average]=0.
        u1_foc_nadir[start:stop, k * n_average:(k + 1) * n_average] = 0.
        u2_foc_nadir[start:stop, k * n_average:(k + 1) * n_average] = 0.


    return num_blanked_lines, blanking_dB_rel_max_amp

def nadir_echo_blanking_thrashold(u_foc_nadir,u1_foc_nadir, u2_foc_nadir, n_average, blanking_dB_over_smeard_signal):
    blanking_value_rel_signal=10**(blanking_dB_over_smeard_signal/20)
    size=u_foc_nadir.shape
    n_times=np.int(np.round(size[1]/n_average))
    num_blanked_lines = np.zeros(n_times)
    for k in range(0, n_times):
        print(k)
        part=u_foc_nadir[8000:16384, k * n_average:(k + 1) * n_average]
        size_part=part.shape
        nadir_shape = np.abs(part).sum(axis=1)/size_part[1]
        signal_level=nadir_shape[150]
        mask = np.where(nadir_shape > signal_level*blanking_value_rel_signal, 1, 0)
        mask_one=np.array(np.where(mask == 1))
        i=mask_one[0,0]
        while mask[i] == 1 and i<size_part[0]-1:
            i = i + 1
        start=mask_one[0,0]
        stop=i
        num_blanked_lines[k]= stop-start
        u_foc_nadir[8000+start:8000+stop,k * n_average:(k + 1) * n_average]=0.
        u1_foc_nadir[8000+start:8000+stop, k * n_average:(k + 1) * n_average] = 0.
        u2_foc_nadir[8000+start:8000+stop, k * n_average:(k + 1) * n_average] = 0.


    return num_blanked_lines, blanking_dB_over_smeard_signal


