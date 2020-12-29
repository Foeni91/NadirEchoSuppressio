import numpy as np
from Raw_data_generation import Azimuth_Modulation
from Raw_data_generation import Range_Modulation


def raw_data_generation(param, system_flag):
    # %%%%%%%%%%%%%     INTRODUCING THE AZIMUTH MODULATION     %%%%%%%%%%%%%

    Azimuth_Modulation.azimuth_modulation_signal(param)

    # %%%%%%%%%%%%%     INTRODUCING THE CHIRP MODULATION     %%%%%%%%%%%%%

    if system_flag == 'WE-SAR':

        Range_Modulation.range_modulation_WESAR_signal(param, system_flag)

    elif system_flag == 'ConvSAR':

        Range_Modulation.range_modulation_conv_signal(param)

    elif system_flag == 'Both':

        Range_Modulation.range_modulation_WESAR_signal(param, system_flag)
        Range_Modulation.range_modulation_conv_signal(param)


