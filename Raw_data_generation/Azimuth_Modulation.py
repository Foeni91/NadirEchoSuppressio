import numpy as np
from scipy.io.idl import readsav
from Raw_data_generation import Azimuth_Antenna_Pattern
import FFT


def azimuth_modulation_signal(param):

    # %%%%%%%%%%%%%     INTRODUCING THE AZIMUTH MODULATION     %%%%%%%%%%%%%

    # Load the data
    if param.signal == 1:
        data_signal = readsav('Scenes/compl_refl_lake.sav')
    if param.signal == 2:
        data_signal = readsav('Scenes/compl_refl_city.sav')
    if param.signal == 3:
        data_signal = readsav('Scenes/compl_refl_forest.sav')
    if param.signal == 4:
        data_signal = readsav('Scenes/compl_refl_town.sav')

    # Scene
    signal_compl_refl_zp = np.zeros((param.n_az, param.n_rg), dtype=complex)
    signal_compl_refl_zp[int(1 * param.n_az / 4):int(3 * param.n_az / 4),
    int(3 * param.n_rg / 8):int(5 * param.n_rg / 8)] = data_signal.compl_refl       # Use for data 16384x8192 (just comment out or in)
    #int(7 * param.n_rg / 16): int(9 * param.n_rg / 16)] = data_signal.compl_refl   # Use for data 32768x8192 (just comment out or in)

    # Azimuth modulation
    FFT.my_ifft_azimuth(FFT.my_fft_azimuth(signal_compl_refl_zp) *
                        FFT.my_fft_azimuth(Azimuth_Antenna_Pattern.antenna_pattern_signal(param)), savpath='u1')

    np.save('u1_conv', np.load('u1.npy', mmap_mode='r'))

