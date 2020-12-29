import Parameters
from Raw_data_generation import Raw_Data_Generation


if __name__ == '__main__':

    ' %%%%%%%%%%%%%     PARAMETERS     %%%%%%%%%%%%% '

    param = Parameters.Parameters()

    #param.print()

    ' %%%%%%%%%%%%%     RAW DATA GENERATION     %%%%%%%%%%%%% '

    # system_flag: 'WE-SAR' OR 'ConvSAR' OR 'Both'
    # (Which system to simulate)
    system_flag = 'WE-SAR'

    Raw_Data_Generation.raw_data_generation(param, system_flag)

    '%%%%%%%%%% Save parameters in Parameters.txt %%%%%'

    f = open('Parameters.txt', "w")
    f.write('# Constants\n'
            'c0 = ' + str(param.c) + '\t\t\t\t# speed of light in vacuum [m/s]\n'
            'R_E = ' + str(param.R_EARTH) + '\t\t\t\t\t# Earths mean radius [m]\n\n'
            'G =' + str(param.G) + '\t\t\t\t\t# Gravitational constant [m^3 / kg / s^2]\n\n'
             'm_E =' + str(param.M_EARTH) + '\t\t\t\t\t# Earth mass [kg]\n\n'
                                            
            '# System parameters\n\n'

            'lamb = ' + str(param.wavelength) + '\t\t\t\t\t# wavelength [m]\n'
            'h = ' + str(param.height) + '\t\t\t\t\t# orbit height [m]\n'
            'CR = ' + str(param.cr) + '\t\t\t\t\t\t# chirp compression ratio\n'
            'B = ' + str(param.chirp_bandwidth) + '\t\t\t\t\t# chirp bandwidth [Hz]\n'
            'tau = ' + str(param.chirp_duration) + '\t# uncompressed pulse length [s] (Calculated: tau = CR/B)\n'
            'f_s = ' + str(param.rg_f_s) + '\t\t\t\t# range sampling frequency [Hz]\n'
            'PRF = ' + str(param.PRF) + '\t\t\t\t\t\t# pulse repetition frequency [Hz]\n'
            'L = ' + str(param.antenna_length) + '\t\t\t\t\t\t\t# antenna length [m]\n'
            'R0= ' + str(param.R0) + '\t\t\t\t\t# target distance [m] (Calculated: R0 = h + c0/2/PRF)\n\n'

            'n_rg = ' + str(param.n_rg) + '\t\t\t\t\t# number of range samples (Calculated: n_rg = 2^np.ceil(np.log(2 * tau * f_s) / np.log(2))\n'
            'n_az = ' + str(param.n_az) + '\t\t\t\t\t\t# number of azimuth samples\n'

            )
    f.close()