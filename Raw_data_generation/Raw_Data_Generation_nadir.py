import numpy as np
from Raw_data_generation import Azimuth_Modulation_nadir
from Raw_data_generation import Range_Modulation_nadir


def raw_data_generation_nadir(param, system_flag):
    # %%%%%%%%%%%%%     INTRODUCING THE AZIMUTH MODULATION     %%%%%%%%%%%%%

    Azimuth_Modulation_nadir.azimuth_modulation_nadir(param)

    # %%%%%%%%%%%%%     INTRODUCING THE CHIRP MODULATION     %%%%%%%%%%%%%

    if system_flag == 'WE-SAR':

        Range_Modulation_nadir.range_modulation_WESAR_nadir(param, system_flag)

    elif system_flag == 'ConvSAR':

        Range_Modulation_nadir.range_modulation_conv_nadir(param)

    elif system_flag == 'Both':

        Range_Modulation_nadir.range_modulation_WESAR_nadir(param, system_flag)
        Range_Modulation_nadir.range_modulation_conv_nadir(param)


