import numpy as np
import math
import FFT


def chirp_modulation_WESAR_signal(param):

    # %%%%%%%%%%%%%     MATRIX INITIALIZATION     %%%%%%%%%%%%%

    # Filter matrix for chirp modulation
    h_rg_sign_mod_matrix = np.zeros((param.n_az, param.n_rg), dtype=complex)

    # %%%%%%%%%%%%%     INTRODUCING THE CHIRP MODULATION     %%%%%%%%%%%%%

    for k in range(0, param.n_az):

        # If cyclically-shifted chirps on transmit
        if param.waveform_sequence == 1:

            t_i_local_sign = param.t_i[k % len(param.t_i)]

            # Range filter for chirp modulation: WE-SAR (focusing matched to the signal)
            h_rg_sign_mod_matrix[k, :] = np.roll(np.where(abs(param.t_rg - 2 * param.R0 / param.c) <= (
                    param.chirp_duration / 2), 1, 0) * np.exp(
                param.j * math.pi * param.chirp_bandwidth / param.chirp_duration *
                (param.t_rg - 2 * param.R0 / param.c - t_i_local_sign - param.chirp_duration * np.floor(
                    (param.t_rg - 2 * param.R0 / param.c + param.chirp_duration / 2 - t_i_local_sign)
                    / param.chirp_duration)) ** 2), int(param.n_rg / 2))

        # If up- down-chirp alternation
        elif param.waveform_sequence == 2:

            up_down_factor_sign = param.up_down_factorVec[k % len(param.up_down_factorVec)]

            # Range filter for chirp modulation: WE-SAR (focusing matched to the signal)
            h_rg_sign_mod_matrix[k, :] = np.roll(
                np.where(abs(param.t_rg - 2 * param.R0 / param.c) <= (param.chirp_duration / 2), 1, 0) * np.exp(
                    param.j * up_down_factor_sign * math.pi * param.chirp_bandwidth / param.chirp_duration * (
                            param.t_rg - 2 * param.R0 / param.c) ** 2), int(param.n_rg / 2))

    return h_rg_sign_mod_matrix

def chirp_modulation_conv(param):

    # %%%%%%%%%%%%%     MATRIX INITIALIZATION     %%%%%%%%%%%%%

    # Filter matrix for chirp modulation
    h_rg_desided_mod_matrix = np.zeros((param.n_az, param.n_rg), dtype=complex)

    # %%%%%%%%%%%%%     INTRODUCING THE CHIRP MODULATION     %%%%%%%%%%%%%

    # Range filter for chirp modulation: Conventional SAR
    h_rg_desired = np.roll(np.where(abs(param.t_rg - 2 * param.R0 / param.c) <= (param.chirp_duration / 2), 1, 0) *
                           np.exp(param.j * math.pi * param.chirp_bandwidth / param.chirp_duration *
                                  (param.t_rg - 2 * param.R0 / param.c) ** 2), int(param.n_rg / 2))

    for k in range(0, param.n_az):

        # Range filter matrix for chirp modulation: Conventional SAR
        h_rg_desided_mod_matrix[k, :] = h_rg_desired

    return h_rg_desided_mod_matrix
