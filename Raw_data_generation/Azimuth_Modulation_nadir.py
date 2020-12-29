import numpy as np
from scipy.io.idl import readsav
from Raw_data_generation import Azimuth_Antenna_Pattern
import FFT


def azimuth_modulation_nadir(param):

    # %%%%%%%%%%%%%     INTRODUCING THE AZIMUTH MODULATION     %%%%%%%%%%%%%

    signal_compl_refl_zp=np.load('nadir_echo.npy', mmap_mode='r').T  # Load the nadir echos signal (backscatter) file nadir_echo.npy generated with nadir_echo_generation.py

    # Azimuth modulation
    FFT.my_ifft_azimuth(FFT.my_fft_azimuth(signal_compl_refl_zp) *
                        FFT.my_fft_azimuth(Azimuth_Antenna_Pattern.antenna_pattern_signal(param)), savpath='nadir_echo_az')