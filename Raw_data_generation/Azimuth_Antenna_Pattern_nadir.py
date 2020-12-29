import math
import numpy as np


def antenna_pattern_nadir(param):
    # %%%%%%%%%%%%%     MATRIX INITIALIZATION     %%%%%%%%%%%%%

    az_mod_sign_matrix = np.zeros((param.n_az, param.n_rg), dtype=complex)

    # %%%%%%%%%%%%%     AZIMUTH ANTENNA PATTERN & MODULATION     %%%%%%%%%%%%%

    az_mod_sign = np.roll(np.sinc(param.antenna_length / param.wavelength * param.sin_psi_nad) ** 2 * np.exp(
        -param.j * 4 * math.pi / param.wavelength * param.h_nad), int(param.n_az / 2))

    #az_mod_sign = np.roll(np.exp(-param.j * 4 * math.pi / param.wavelength * param.R), int(param.n_az / 2))

    for k in range(0, param.n_rg):
        az_mod_sign_matrix[:, k] = az_mod_sign

    return az_mod_sign_matrix