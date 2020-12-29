"""
# author Fabian Glatz
##################################################################
# Description:
# This code calculates the total error (Error due to the smeared nadir echo + corruption of the useful signal due to the blanking)
# as a function of the blanked lines. This function has a minimum. When the errors are equally weighted, the number of blanked lines
# at this minimum leads to the best achievable image quality.
# The values of the function are stored in a file.
# The simulation in done for a waveform-encoded SAR using shifted chirps.
##################################################################
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal.windows import general_hamming
import Parameters
import os
from datetime import datetime
import nadir_echo_blanking

# DEMOSTRATION OF NADIR ECHO SUPPRESSION WITH A REAL TERRASAR-X SCENE AND A REALISTIC NADIR ECHO USING SHIFTED CHIRPS

param = Parameters.Parameters()         # Load constants and parameters from parameters class

# Constants

c0 = param.c                            # speed of light in vacuum [m/s]
R_E = param.R_EARTH                     # Earth's mean radius [m]

# System parameters

lamb = param.wavelength                     # wavelength [m]
h = param.height                            # orbit height [m]
CR = param.cr                               # chirp compression ratio
B = param.chirp_bandwidth                   # chirp bandwidth [Hz]
tau = param.chirp_duration                  # uncompressed pulse length [s]
f_s = param.rg_f_s                          # range sampling frequency [Hz]
PRF = param.PRF                             # pulse repetition frequency [Hz]
L = param.antenna_length                    # antenna length [m]
R0=param.R0                                 # target distance [m]

# Processing parameters

hamming_rg = param.hamming_rg                       # Hamming window in range (0 -> no, 1 -> yes)
hamming_rg_alpha = param.hamming_rg_alpha           # alpha coefficient for the Hamming window in range
PBW = param.PBW                                     # processed Doppler bandwidth [Hz]
hamming_az = param.hamming_az                       # Hamming window in azimuth (0 -> no, 1 -> yes)
hamming_az_alpha = param.hamming_az_alpha           # alpha coefficient for the Hamming window in azimuth
pattern_az = param.pattern_az                       # azimuth antenna pattern compensation (0 -> no, 1 -> yes)
plots = param.plots                                 # print the plots (0 -> no, 1 -> yes)
data = param.data                                   # print the data (0 -> no, 1 -> yes)

# Code (i.e., sequence of STSO waveforms) to be used on transmit
# Each waveform is characterised by a t_i
# All t_i are collected in a vector

t_i=param.t_i

j = param.j
G = param.G
m_E = param.M_EARTH
v_s = param.v_s
v_g = param.v_g                                                             # ground velocity [m/s]
v_r = param.v_r                                                             # effective velocity [m/s]
# look_angle_main = slant_range_to_look_angle (R0, R_E, h)
n_az = param.n_az                                                           # number of azimuth samples  - in the letter it was n_az = 2.**14.
az_spacing = v_g / PRF                                                      # azimuth spacing [m]
t_az = param.t_az                                                           # azimuth time [s]
n_rg = param.n_rg                                                           # number of range samples
rg_spacing = c0 / 2 / f_s                                                   # slant range spacing [m]
t_rg = param.t_rg                                                           # range time [s]
R = param.R                                                                 # range [m]
R_nad = h + v_r ** 2. * t_az** 2. / (2. * h) + c0/PRF/2.                    # range [m]
h_nad = h + v_r ** 2. * t_az** 2. / (2. * h)
sin_psi = param.sin_psi                                                     # sinus of the azimuth angle
sin_psi_nad = (v_g * t_az / h) / np.sqrt(1. + (v_g * t_az / h) ** 2.)       # sinus of the azimuth angle
psi = np.arcsin(sin_psi)
psi_nad = np.arcsin(sin_psi_nad)

signal_equiv_peak=1.0   # Selection of the averaged maximum backscatter sigma_max_0 of the nadir echo.
                        # This value always references the value of the backscatter that would result if a conventional SAR system were used.


# ********************************************************
# Generate Folder and save Parameters
# ********************************************************
# datetime object containing current date and time
now = datetime.now()
# dd/mm/YY H:M:S
dt_string = now.strftime("%Y%m%d_%H%M%S")
folder_name = 'Shifted_Chirp_no_RCM_CR_' + str(CR) + '_n_az_2_' + str(int(np.log2(n_az))) + '_' + dt_string
os.mkdir(folder_name)
f = open(folder_name + '\Parameters.txt', "w")
f.write('# Constants\n'
        'c0 = ' + str(c0) + '\t\t\t\t# speed of light in vacuum [m/s]\n'
        'R_E = ' + str(R_E) + '\t\t\t\t\t# Earths mean radius [m]\n'
        'G =' + str(G) + '\t\t\t\t\t# Gravitational constant [m^3 / kg / s^2]\n'
        'm_E =' + str(m_E) + '\t\t\t\t\t# Earth mass [kg]\n\n'
        '# System parameters\n\n'
        'lamb = ' + str(lamb) + '\t\t\t\t\t# wavelength [m]\n'
        'h = ' + str(h) + '\t\t\t\t\t# orbit height [m]\n'
        'CR = ' + str(CR) + '\t\t\t\t\t\t# chirp compression ratio\n'
        'B = ' + str(B) + '\t\t\t\t\t# chirp bandwidth [Hz]\n'
        'tau = ' + str(tau) + '\t# uncompressed pulse length [s] (Calculated: tau = CR/B)\n'
        'f_s = ' + str(f_s) + '\t\t\t\t# range sampling frequency [Hz]\n'
        'PRF = ' + str(PRF) + '\t\t\t\t\t\t# pulse repetition frequency [Hz]\n'
        'L = ' + str(L) + '\t\t\t\t\t\t\t# antenna length [m]\n'
        'R0= ' + str(R0) + '\t\t\t\t\t# target distance [m] (Calculated: R0 = h + c0/2/PRF)\n\n'
        '# Processing parameters\n\n'
        'hamming_rg = ' + str(hamming_rg) + '\t\t\t\t# Hamming window in range (0 -> no, 1 -> yes)\n'
        'hamming_rg_alpha = ' + str(hamming_rg_alpha) + '\t\t\t# alpha coefficient for the Hamming window in range\n'
        'PBW = ' + str(PBW) + '\t\t\t\t\t\t# processed Doppler bandwidth [Hz]\n'
        'hamming_az = ' + str(hamming_az) + '\t\t\t\t# Hamming window in azimuth (0 -> no, 1 -> yes)\n'
        'hamming_az_alpha = ' + str(hamming_az_alpha) + '\t\t\t# alpha coefficient for the Hamming window in azimuth\n'
        'pattern_az = ' + str(pattern_az) + '\t\t\t\t# azimuth antenna pattern compensation (0 -> no, 1 -> yes)\n'
        'plots = ' + str(plots) + '\t\t\t\t\t\t# print the plots (0 -> no, 1 -> yes)\n'
        'data = ' + str(data) + '\t\t\t\t\t\t# print the data (0 -> no, 1 -> yes)\n\n'

        'n_rg = ' + str(n_rg) + '\t\t\t\t\t# number of range samples (Calculated: n_rg = 2^np.ceil(np.log(2 * tau * f_s) / np.log(2))\n'
        'n_az = ' + str(n_az) + '\t\t\t\t\t\t# number of azimuth samples\n'

        'Nadir echo peak = ' + str(signal_equiv_peak) + '\t\t\t\t\t\t# Peak of the nadir echo in the focused image\n'

        )
f.close()

u1=np.load('Raw_data_generation/u1_raw_shifted_chirp_rg_az.npy')            # load raw data of the scene
u1=u1.T*0.047584658951353415                                                # The calibration values must be chosen so that the processed data of the scenes correspond to the real ones.
                                                                            # Real values see Nertjanas investigation / Calibration values for this code see Calibraion.txt in the project folder


# *****

# Filtering adapted to the signal

# *****

# RANGE COMPRESSION

if (hamming_rg):

  n_rg_eff = 2. * np.round((n_rg * B / f_s)/2.)

  h_hamming_rg_F = np.roll(np.concatenate((np.repeat(0., (n_rg - n_rg_eff) / 2.), general_hamming(int(n_rg_eff),hamming_rg_alpha), np.repeat(0., (n_rg - n_rg_eff) / 2.)),axis=0), int(n_rg / 2))

# Hamming window in range

h_rg_desired = np.roll(np.less_equal(np.abs(t_rg-2*R0/c0), tau/2) *np.exp(- j* np.pi * B / tau * (t_rg - 2. * R0 / c0) ** 2.), int(n_rg / 2.))
# impulse response of the range compression filter

temp_var=np.fft.fft(h_rg_desired)*np.conj(np.fft.fft(h_rg_desired))

for k in range(0, n_az):
    t_i_local = t_i[(k) % len(t_i)]

    # h_rg = np.roll((np.abs(t_rg - 2. * R0 / c0) <= (tau / 2.)) * np.exp(j * np.pi * B / tau * (t_rg - 2. * R0 / c0 - t_i_local - tau * np.floor((t_rg - 2. * R0 / c0 + tau / 2 - t_i_local) / tau)) ** 2.), int(n_rg / 2.))

    h_rg = np.roll(np.less_equal(np.abs(t_rg - 2 * R0 / c0), tau / 2) * np.exp(j * np.pi * B / tau * (
                t_rg - 2. * R0 / c0 - t_i_local - tau * np.floor(
            (t_rg - 2. * R0 / c0 + tau / 2 - t_i_local) / tau)) ** 2.), int(n_rg / 2.))

    # h_rg_F_ideal=temp_var/np.fft.fft(h_rg)
    h_rg_F_ideal = temp_var / np.fft.fft(h_rg)

    if (hamming_rg):
        h_rg_F_ideal = h_rg_F_ideal * h_hamming_rg_F

    u1[:, k] = np.fft.ifft(np.fft.fft(u1[:, k]) * h_rg_F_ideal)     # range compression (performed in the frequency domain)








f_rg = np.roll((np.arange(n_rg) - n_rg / 2.) * f_s / n_rg, int(n_rg / 2.)) # range frequency vector [Hz]
f_az = np.roll((np.arange(n_az) - n_az / 2.) * PRF / n_az, int(n_az / 2.)) # azimuth frequency vector [Hz]
f0 = c0 / lamb    # radar carrier frequency [Hz]

h_az = np.exp(1j * 4. * np.pi / lamb * R)
t_p = float(R0 * lamb / L / v_g)  # observation time [s]
h_az[np.where(abs(t_az) >= (t_p / 2))] = 0
h_az = np.roll(h_az, int(n_az / 2))
h_az_F = np.fft.fft(h_az)
for k in range(0, n_rg):
    u1[k, :] = np.fft.fft(u1[k, :]) * h_az_F


n_az_eff = 2. * round((n_az*PBW/PRF)/2.)

PBW_window_F =np.roll(np.concatenate((np.repeat(0., (n_az - n_az_eff) / 2.), np.repeat (1., n_az_eff), np.repeat(0.,(n_az - n_az_eff)/2.)),axis=0), int(n_az / 2)) # azimuth window to select the processed Doppler bandwidth

az_window_F = PBW_window_F

if (hamming_az):
    h_hamming_az_F = np.roll(np.concatenate((np.repeat(0., (n_az - n_az_eff) / 2.), general_hamming(int(n_az_eff), hamming_az_alpha), np.repeat(0.,(n_az - n_az_eff)/2.)),axis=0), int(n_az/2))
    az_window_F = az_window_F * h_hamming_az_F
    # Hamming window in azimuth

if (pattern_az):
    h_az_pattern_comp_F = (abs(f_az) <= (PBW / 2.)) / np.sinc(L / v_g / 2. * f_az) ** 2.
    az_window_F = az_window_F * h_az_pattern_comp_F
    # compensaion of the two-way aziuth antnna pattern

for k in range(0, n_rg):
    u1[k, :] = np.fft.ifft(u1[k, :] * az_window_F)



u1_single=u1.copy()
np.save(folder_name+'/u1_single', u1_single)









u1=np.load('Raw_data_generation/u1_raw_shifted_chirp_rg_az.npy')            # load raw data of the scene
u1=u1.T*0.047584658951353415                                                # The calibration values must be chosen so that the processed data of the scenes correspond to the real ones.
                                                                            # Real values see Nertjanas investigation / Calibration values for this code see Calibraion.txt in the project folder
u2=np.load('Raw_data_generation/nadir_echo_raw_shifted_chirp_rg_az.npy')    # load raw data of the nadir echo
u2=u2.T/7.801635488463321*234.18157380078924*signal_equiv_peak              # These calibration values lead to a maximum backscatter of the nadir echo of 1.
                                                                            # Therefore, the maximum backsactter can be freely selected via signal_equiv_peak.


data_values=open(folder_name+'\data.txt', 'w')
print( 'Raw data')
data_values.write('Raw data\n\n')






# *****
# Filtering adapted to the nadir echo
# *****

# RANGE COMPRESSION

if (hamming_rg):

  n_rg_eff = 2. * np.round((n_rg * B / f_s)/2.)

  h_hamming_rg_F = np.roll(np.concatenate((np.repeat(0., (n_rg - n_rg_eff) / 2.), general_hamming(int(n_rg_eff),hamming_rg_alpha), np.repeat(0., (n_rg - n_rg_eff) / 2.)),axis=0), int(n_rg / 2))
# Hamming window in range

h_rg_desired = np.roll(np.less_equal(np.abs(t_rg-2*R0/c0), tau/2) *np.exp(- j* np.pi * B / tau * (t_rg - 2. * R0 / c0) ** 2.), int(n_rg / 2.))
# # impulse response of the range compression filter

temp_var=np.fft.fft(h_rg_desired)*np.conj(np.fft.fft(h_rg_desired))

for k in range(0,n_az):

    t_i_local=t_i[(k+1) % len(t_i)]

    #h_rg = np.roll((np.abs(t_rg - 2. * R0 / c0) <= (tau / 2.)) * np.exp(j * np.pi * B / tau * (t_rg - 2. * R0 / c0 - t_i_local - tau * np.floor((t_rg - 2. * R0 / c0 + tau / 2 - t_i_local) / tau)) ** 2.), int(n_rg / 2.))
    h_rg = np.roll(np.less_equal(np.abs(t_rg - 2 * R0 / c0), tau / 2) * np.exp(j * np.pi * B / tau * (t_rg - 2. * R0 / c0 - t_i_local - tau * np.floor((t_rg - 2. * R0 / c0 + tau / 2 - t_i_local) / tau)) ** 2.), int(n_rg / 2.))

    h_rg_F_ideal = temp_var / np.fft.fft(h_rg)

    if (hamming_rg):
        h_rg_F_ideal = h_rg_F_ideal * h_hamming_rg_F


    u1[:, k] = np.fft.ifft(np.fft.fft(u1[:, k]) * h_rg_F_ideal)  # range compression (performed in the frequency domain)
    u2[:, k] = np.fft.ifft(np.fft.fft(u2[:, k]) * h_rg_F_ideal)  # range compression (performed in the frequency domain)



print( 'Range-compressed data')
data_values.write('Range compressed data to the nadir echo\n\n')


# RANGE CELL MIGRATION (RCM) CORRECTION AND AZIMUTH COMPRESSION



f_rg = np.roll((np.arange(n_rg) - n_rg / 2.) * f_s / n_rg, int(n_rg / 2.)) # range frequency vector [Hz]
f_az = np.roll((np.arange(n_az) - n_az / 2.) * PRF / n_az, int(n_az / 2.)) # azimuth frequency vector [Hz]
f0 = c0 / lamb    # radar carrier frequency [Hz]


h_az = np.exp(1j * 4. * np.pi / lamb * h_nad)
t_p = float(R0 * lamb / L / v_g)  # observation time [s]
h_az[np.where(abs(t_az) >= (t_p / 2))] = 0
h_az = np.roll(h_az, int(n_az / 2))
h_az_F = np.fft.fft(h_az)
for k in range(0, n_rg):
 ############# but you have to make a cycle for all ranges
    u1[k, :] = np.fft.ifft(np.fft.fft(u1[k, :]) * h_az_F)
    u2[k, :] = np.fft.ifft(np.fft.fft(u2[k, :]) * h_az_F)

n_az_eff = 2. * round((n_az*PBW/PRF)/2.)

PBW_window_F =np.roll(np.concatenate((np.repeat(0., (n_az - n_az_eff) / 2.), np.repeat (1., n_az_eff), np.repeat(0.,(n_az - n_az_eff)/2.)),axis=0), int(n_az / 2)) # azimuth window to select the processed Doppler bandwidth

az_window_F = PBW_window_F




# u = np.fft.ifft(u, axis=1)
# u1 = np.fft.ifft(u1, axis=1)
# u2 = np.fft.ifft(u2, axis=1)




print( 'Focused data')
data_values.write('Focused data to the nadir echo\n\n')



# Save u1 and u2 to new variable

u1_focused_nadir=u1.copy()
u2_focused_nadir=u2.copy()

i=np.arange(0,100,10)                   # Definition how much lines shall be blanked per iteration and how much in total
error=np.zeros([4,len(i)])              # generate a matrix for the data
x=0
for i in i:
    # *****
    # Removal of the nadir echo
    # *****


    data_values.write('u1[8191+4:8191+4' +str(i)+', :] = 0.\n')
    data_values.write('u2[8191+4:8191+4'+str(i)+', :] = 0.\n\n\n')


    u1[8191+4:8191+i+4, :] = 0.
    u2[8191+4:8191+i+4, :] = 0.



    print( 'Focused data after nadir removal')
    data_values.write('Focused data to the nadir echo after nadir echo removal\n\n')


    h_az = np.exp(1j * 4. * np.pi / lamb * h_nad)
    t_p = float(R0 * lamb / L / v_g)  # observation time [s]
    h_az[np.where(abs(t_az) >= (t_p / 2))] = 0
    h_az = np.roll(h_az, int(n_az / 2))
    h_az_F = np.fft.fft(h_az)
    for k in range(0, n_rg):
     ############# but you have to make a cycle for all ranges
        u1[k, :] = np.fft.ifft(np.fft.fft(u1[k, :]) / h_az_F)
        u2[k, :] = np.fft.ifft(np.fft.fft(u2[k, :]) / h_az_F)


    # End test to speed up #


    print( 'Range-compressed data')
    data_values.write('Data after de-compression\n\n')



    # RANGE "DECOMPRESSION"
    for k in range(0, n_az):

        t_i_local = t_i[(k + 1) % len(t_i)]

        # h_rg = np.roll((np.abs(t_rg - 2 * R0 / c0) <= (tau / 2)) * np.exp(j * np.pi * B / tau * (t_rg - 2. * R0 / c0 - t_i_local - tau * np.floor((t_rg - 2. * R0 / c0 + tau / 2 - t_i_local) / tau)) ** 2.), int(n_rg / 2.))
        h_rg = np.roll(np.less_equal(np.abs(t_rg - 2 * R0 / c0), tau / 2) * np.exp(j * np.pi * B / tau * (
                    t_rg - 2. * R0 / c0 - t_i_local - tau * np.floor(
                (t_rg - 2. * R0 / c0 + tau / 2 - t_i_local) / tau)) ** 2.), int(n_rg / 2.))

        h_rg_F_ideal = temp_var / np.fft.fft(h_rg)

        if (hamming_rg):
            h_rg_F_ideal = h_rg_F_ideal * h_hamming_rg_F

        # h_rg_F_ideal_mask = np.where(h_rg_F_ideal == 0.,1., 0)
        n_rg_eff = 2. * round((n_rg * B / f_s) / 2.)

        h_rg_F_ideal_mask = np.roll(np.concatenate(
            (np.repeat(1., (n_rg - n_rg_eff) / 2.), np.repeat(0., n_rg_eff), np.repeat(1., (n_rg - n_rg_eff) / 2.)),
            axis=0), int(n_rg / 2))  # azimuth window to select the processed Doppler bandwidth


        u1[:, k] = np.fft.ifft(np.fft.fft(u1[:, k]) * 1. / (h_rg_F_ideal + h_rg_F_ideal_mask) * (
                    1. - h_rg_F_ideal_mask))  # range compression (performed in the frequency domain)
        u2[:, k] = np.fft.ifft(np.fft.fft(u2[:, k]) * 1. / (h_rg_F_ideal + h_rg_F_ideal_mask) * (
                    1. - h_rg_F_ideal_mask))  # range compression (performed in the frequency domain)




    print( 'Raw data')
    data_values.write('Raw data after defocusing\n\n')

    # *****
    # Filtering adapted to the signal
    # *****

    # RANGE COMPRESSION
    for k in range(0, n_az):

        t_i_local = t_i[k % len(t_i)]

        # h_rg = np.roll((np.abs(t_rg - 2 * R0 / c0) <= (tau / 2)) * np.exp(j * np.pi * B / tau * (t_rg - 2. * R0 / c0 - t_i_local - tau * np.floor((t_rg - 2. * R0 / c0 + tau / 2 - t_i_local) / tau)) ** 2.), int(n_rg / 2.))
        h_rg = np.roll(np.less_equal(np.abs(t_rg - 2 * R0 / c0), tau / 2) * np.exp(j * np.pi * B / tau * (
                    t_rg - 2. * R0 / c0 - t_i_local - tau * np.floor(
                (t_rg - 2. * R0 / c0 + tau / 2 - t_i_local) / tau)) ** 2.), int(n_rg / 2.))

        h_rg_F_ideal = temp_var / np.fft.fft(h_rg)

        if (hamming_rg):
            h_rg_F_ideal = h_rg_F_ideal * h_hamming_rg_F


        u1[:, k] = np.fft.ifft(
            np.fft.fft(u1[:, k]) * h_rg_F_ideal)  # range compression (performed in the frequency domain)
        u2[:, k] = np.fft.ifft(
            np.fft.fft(u2[:, k]) * h_rg_F_ideal)  # range compression (performed in the frequency domain)


    print('Range-compressed data')
    data_values.write('Range compressed data to the signal after nadir echo removal \n\n')


    f_rg = np.roll((np.arange(n_rg) - n_rg / 2.) * f_s / n_rg, int(n_rg / 2.)) # range frequency vector [Hz]
    f_az = np.roll((np.arange(n_az) - n_az / 2.) * PRF / n_az, int(n_az / 2.)) # azimuth frequency vector [Hz]
    f0 = c0 / lamb    # radar carrier frequency [Hz]


    h_az = np.exp(1j * 4. * np.pi / lamb * R)
    t_p = float(R0 * lamb / L / v_g)  # observation time [s]
    h_az[np.where(abs(t_az) >= (t_p / 2))] = 0
    h_az = np.roll(h_az, int(n_az / 2))
    h_az_F = np.fft.fft(h_az)
    for k in range(0, n_rg):
     ############# but you have to make a cycle for all ranges
        u1[k, :] = np.fft.fft(u1[k, :]) * h_az_F
        u2[k, :] = np.fft.fft(u2[k, :]) * h_az_F

    n_az_eff = 2. * round((n_az*PBW/PRF)/2.)

    PBW_window_F =np.roll(np.concatenate((np.repeat(0., (n_az - n_az_eff) / 2.), np.repeat (1., n_az_eff), np.repeat(0.,(n_az - n_az_eff)/2.)),axis=0), int(n_az / 2)) # azimuth window to select the processed Doppler bandwidth

    az_window_F = PBW_window_F

    if (hamming_az):
        h_hamming_az_F = np.roll(np.concatenate((np.repeat(0., (n_az - n_az_eff) / 2.), general_hamming(int(n_az_eff), hamming_az_alpha), np.repeat(0.,(n_az - n_az_eff)/2.)),axis=0), int(n_az/2))
        az_window_F = az_window_F * h_hamming_az_F
        # Hamming window in azimuth

    if (pattern_az):
        h_az_pattern_comp_F = (abs(f_az) <= (PBW / 2.)) / np.sinc(L / v_g / 2. * f_az) ** 2.
        az_window_F = az_window_F * h_az_pattern_comp_F
        # compensaion of the two-way aziuth antnna pattern

    for k in range(0, n_rg):
        u1[k, :] = np.fft.ifft(u1[k, :] * az_window_F)
        u2[k, :] = np.fft.ifft(u2[k, :] * az_window_F)

    print('Focused data')
    data_values.write('Focused data to the signal after nadir echo removal \n\n')

    error_nadir_echo = (np.abs(u2[6144:10240, 2048:6144]) ** 2).sum()                                       # calculate the error causes be the nadir echo
    error_signal = (np.abs(u1_single[6144:10240, 2048:6144] - u1[6144:10240, 2048:6144]) ** 2).sum()        # calculate the image corruption due to the blanking
    error_tot = (np.abs(u2[6144:10240, 2048:6144]) ** 2).sum() + (
                np.abs(u1_single[6144:10240, 2048:6144] - u1[6144:10240, 2048:6144]) ** 2).sum()            # calculate the total error

    error[0,x]=i                                    # Write the values in a data table
    error[1, x] = error_nadir_echo                  # Write the values in a data table
    error[2, x] = error_signal                      # Write the values in a data table
    error[3, x] = error_tot                         # Write the values in a data table

    x=x+1
    u1=u1_focused_nadir.copy()
    u2=u2_focused_nadir.copy()


np.savetxt(folder_name+'\error.cvs',error,delimiter=',')
data_values.close()
