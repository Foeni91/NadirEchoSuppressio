import numpy as np
import FFT
from Raw_data_generation import Filters


def range_modulation_WESAR_nadir(param, system_flag):

    # %%%%%%%%%%%%%     WE-SAR: LOAD AZIMUTH MODULATED DATA     %%%%%%%%%%%%%

    u1 = np.load('nadir_echo_az.npy', mmap_mode='r')

    # %%%%%%%%%%%%%     WE-SAR: INTRODUCING THE CHIRP MODULATION     %%%%%%%%%%%%%

    # Range modulation (performed in the frequency domain): WE-SAR
    if param.waveform_sequence == 1:
        filter_mod = Filters.chirp_modulation_WESAR_signal(param)
        filter_mod = np.roll(filter_mod, -1, axis=0)
        FFT.my_ifft_range(FFT.my_fft_range(u1) * FFT.my_fft_range(filter_mod), savpath='nadir_echo_raw_shifted_chirp_rg_az')
    elif param.waveform_sequence == 2:
        filter_mod = Filters.chirp_modulation_WESAR_signal(param)
        filter_mod=np.roll(filter_mod, -1, axis=0)
        FFT.my_ifft_range(FFT.my_fft_range(u1) * FFT.my_fft_range(filter_mod), savpath='nadir_echo_raw_up_down_chirp_rg_az')


def range_modulation_conv_nadir(param):

    # %%%%%%%%%%%%%     CONVENTIONAL SAR: LOAD AZIMUTH MODULATED DATA     %%%%%%%%%%%%%

    u1_conv = np.load('nadir_echo_az.npy', mmap_mode='r')

    # %%%%%%%%%%%%%     INTRODUCING THE CHIRP MODULATION     %%%%%%%%%%%%%

    # Range modulation (performed in the frequency domain): Conventional SAR
    filter_mod = Filters.chirp_modulation_conv(param)
    FFT.my_ifft_range(FFT.my_fft_range(u1_conv) * FFT.my_fft_range(filter_mod), savpath='nadir_echo_raw_conv_rg_az')