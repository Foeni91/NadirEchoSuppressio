"""
# author Fabian Glatz
##################################################################
# Description:
# This code simulates a waveform-encoded SAR system (no range cell migration) using shifted chirps.
# A real scene is overlaid with a previously generated real nadir echo signal.
# The SAR raw data for the scene and the nadir echo have to be generated in advance using the code main_raw_data_generation.py
# (for the scene) and main_raw_data_generation_nadir.py (for the nadir echo)
# Code for the simulation of a data set N_rg x N_az = 16384 x 8192 (This is just relevant for the blanking of the nadir echo and imaging reasons.)
##################################################################
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal.windows import general_hamming
import Parameters
import os
from datetime import datetime


# DEMOSTRATION OF NADIR ECHO SUPPRESSION WITH A REAL TERRASAR-X SCENE AND A REALISTIC NADIR ECHO USING SHIFTED CHIRPS

param = Parameters.Parameters()         # Load constants and parameters from parameters class

# Constants

c0 = param.c                                                # speed of light in vacuum [m/s]
R_E = param.R_EARTH                                         # Earth's mean radius [m]

# System parameters

lamb = param.wavelength                                     # wavelength [m]
h = param.height                                            # orbit height [m]
CR = param.cr                                               # chirp compression ratio
B = param.chirp_bandwidth                                   # chirp bandwidth [Hz]
tau = param.chirp_duration                                  # uncompressed pulse length [s]
f_s = param.rg_f_s                                          # range sampling frequency [Hz]
PRF = param.PRF                                             # pulse repetition frequency [Hz]
L = param.antenna_length                                    # antenna length [m]
R0=param.R0                                                 # target distance [m]

# Processing parameters

hamming_rg =param.hamming_rg                                # Hamming window in range (0 -> no, 1 -> yes)
hamming_rg_alpha =param.hamming_rg_alpha                    # alpha coefficient for the Hamming window in range
PBW = param.PBW                                             # processed Doppler bandwidth [Hz]
hamming_az = param.hamming_az                               # Hamming window in azimuth (0 -> no, 1 -> yes)
hamming_az_alpha = param.hamming_az_alpha                   # alpha coefficient for the Hamming window in azimuth
pattern_az = param.pattern_az                               # azimuth antenna pattern compensation (0 -> no, 1 -> yes)
plots = param.plots                                         # print the plots (0 -> no, 1 -> yes)
data = param.data                                           # print the data (0 -> no, 1 -> yes)

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

# ********************************************************
# ONLY WAVEFORM VARIATION - NO MULTI-FOCUS POST-PROCESSING
# ********************************************************

u1=np.load('Raw_data_generation/u1_raw_shifted_chirp_rg_az.npy')            # load raw data of the scene
u1=u1.T*0.047584658951353415                                                # The calibration values must be chosen so that the processed data of the scenes correspond to the real ones.
                                                                            # Real values see Nertjanas investigation / Calibration values for this code see Calibraion.txt in the project folder
u2=np.load('Raw_data_generation/nadir_echo_raw_shifted_chirp_rg_az.npy')    # load raw data of the nadir echo
u2=u2.T/7.801635488463321*234.18157380078924*signal_equiv_peak              # These calibration values lead to a maximum backscatter of the nadir echo of 1.
                                                                            # Therefore, the maximum backsactter can be freely selected via signal_equiv_peak.

u=u1+u2                                                                     # Superposition of scene and nadir echo.

# np.save(folder_name+'/u_single_raw', u)
# np.save(folder_name+'/u1_single_raw', u1)
# np.save(folder_name+'/u2_single_raw', u2)

# Plot the data if desired
if (plots):
    # Raw data display
    plt.figure()
    plt.imshow(20*np.log10(np.abs(u)/float(np.abs(u).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Raw_data_u')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Raw_data_u')
    plt.show()

    plt.figure()
    plt.imshow(20*np.log10(np.abs(u1)/float(np.abs(u1).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Raw_data_u')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Raw_data_u1')
    plt.show()

    plt.figure()
    plt.imshow(20*np.log10(np.abs(u2)/float(np.abs(u2).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Raw_data_u')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Raw_data_u2')
    plt.show()

data_values=open(folder_name+'\data.txt', 'w')
data_values.write('Raw data\n\n')
print( 'Raw data')

# Calculate values if desired
if (data):
    data_values.write( 'Peak u1 = '+ str(20.*np.log10(abs(u1).max()))+ ' dB\n')
    data_values.write('Peak u2 = ' + str(20. * np.log10(abs(u2).max())) + ' dB\n')
    data_values.write('Peak u= ' + str(20. * np.log10(abs(u).max())) + ' dB\n\n')

    data_values.write( 'Energy u1= '+ str(10.*np.log10((abs(u1)**2.).sum()))+ ' dB\n')
    data_values.write('Energy u2= ' + str(10. * np.log10((abs(u2) ** 2.).sum())) + ' dB\n')
    data_values.write('Energy u= ' + str(10. * np.log10((abs(u) ** 2.).sum())) + ' dB\n\n\n')

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

for k in range(0,n_az):

    t_i_local=t_i[(k) % len(t_i)]

    h_rg = np.roll(np.less_equal(np.abs(t_rg-2*R0/c0), tau/2) *np.exp(j * np.pi * B / tau * (t_rg - 2. * R0 / c0 - t_i_local - tau * np.floor((t_rg - 2. * R0 / c0 + tau / 2 - t_i_local) / tau)) ** 2.), int(n_rg / 2.))

    h_rg_F_ideal = temp_var / np.fft.fft(h_rg)


    if (hamming_rg):
        h_rg_F_ideal = h_rg_F_ideal * h_hamming_rg_F


    u[:, k] = np.fft.ifft(np.fft.fft(u[:, k]) * h_rg_F_ideal)       # range compression (performed in the frequency domain)
    u1[:, k] = np.fft.ifft(np.fft.fft(u1[:, k]) * h_rg_F_ideal)     # range compression (performed in the frequency domain)
    u2[:, k] = np.fft.ifft(np.fft.fft(u2[:, k]) * h_rg_F_ideal)     # range compression (performed in the frequency domain)


# Plot the data if desired
if (plots):
    plt.figure()
    plt.imshow(20*np.log10(np.abs(u)/float(np.abs(u).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Range compressed data u')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Range_compressed_data_u_to_the_signal_without_nadir_removal')
    plt.show()

    plt.figure()
    plt.imshow(20*np.log10(np.abs(u1)/float(np.abs(u1).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Range compressed data u1')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Range_compressed_data_u1_to_the_signal_without_nadir_removal')
    plt.show()

    plt.figure()
    plt.imshow(20*np.log10(np.abs(u2)/float(np.abs(u2).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Range compressed data u2')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Range_compressed_data_u2_to_the_signal_without_nadir_removal')
    plt.show()

# np.save(folder_name+'/u_single_rangecompressed', u)
# np.save(folder_name+'/u1_single_rangecompressed', u1)
# np.save(folder_name+'/u2_single_rangecompressed', u2)

print( 'Range-compressed data')
data_values.write('Range_compressed data\n\n')

# Calculate values if desired
if (data):
    data_values.write( 'Peak u1 = '+ str(20.*np.log10(abs(u1).max()))+ ' dB\n')
    data_values.write('Peak u2 = ' + str(20. * np.log10(abs(u2).max())) + ' dB\n')
    data_values.write('Peak u= ' + str(20. * np.log10(abs(u).max())) + ' dB\n\n')

    data_values.write( 'Energy u1= '+ str(10.*np.log10((abs(u1)**2.).sum()))+ ' dB\n')
    data_values.write('Energy u2= ' + str(10. * np.log10((abs(u2) ** 2.).sum())) + ' dB\n')
    data_values.write('Energy u= ' + str(10. * np.log10((abs(u) ** 2.).sum())) + ' dB\n\n\n')

#
# AZIMUTH COMPRESSION

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
    u[k, :] = np.fft.fft(u[k, :]) * h_az_F
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
    u[k, :] = np.fft.ifft(u[k, :] * az_window_F)
    u1[k, :] = np.fft.ifft(u1[k, :] * az_window_F)
    u2[k, :] = np.fft.ifft(u2[k, :] * az_window_F)


# Plot the data if desired
if (plots):
    plt.figure()
    plt.imshow(20*np.log10(np.abs(u)/float(np.abs(u).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Focused data u to the signal')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Focused_data_u_to_the_signal_without_dual_focusing')
    plt.show()

    plt.figure()
    plt.imshow(20*np.log10(np.abs(u1)/float(np.abs(u1).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Focused data u1 to the signal')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Focused_data_u1_to_the_signal_without_dual_focusing')
    plt.show()

    plt.figure()
    plt.imshow(20*np.log10(np.abs(u2)/float(np.abs(u2).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Focused data u2 to the signal')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Focused_data_u2_to_the_signal_without_dual_focusing')
    plt.show()


print( 'Focused data')
data_values.write('Focused data to the signal\n\n')


# Calculate values if desired
if (data):
    data_values.write( 'Peak u1 = '+ str(20.*np.log10(abs(u1).max()))+ ' dB\n')
    data_values.write('Peak u2 = ' + str(20. * np.log10(abs(u2).max())) + ' dB\n')
    data_values.write('Peak u= ' + str(20. * np.log10(abs(u).max())) + ' dB\n\n')

    data_values.write( 'Energy u1= '+ str(10.*np.log10((abs(u1)**2.).sum()))+ ' dB\n')
    data_values.write('Energy u2= ' + str(10. * np.log10((abs(u2) ** 2.).sum())) + ' dB\n')
    data_values.write('Energy u= ' + str(10. * np.log10((abs(u) ** 2.).sum())) + ' dB\n\n\n')


# Save processed data after single focusing u, u1, u2
np.save(folder_name+'/u_single_focused', u)
np.save(folder_name+'/u1_single_focused', u1)
np.save(folder_name+'/u2_single_focused', u2)

# **************************************************
# WAVEFORM VARIATION AND MULTI-FOCUS POST-PROCESSING
# **************************************************

u1=np.load('Raw_data_generation/u1_raw_shifted_chirp_rg_az.npy')            # load raw data of the scene
u1=u1.T*0.047584658951353415                                                # The calibration values must be chosen so that the processed data of the scenes correspond to the real ones.
                                                                            # Real values see Nertjanas investigation / Calibration values for this code see Calibraion.txt in the project folder
u2=np.load('Raw_data_generation/nadir_echo_raw_shifted_chirp_rg_az.npy')    # load raw data of the nadir echo
u2=u2.T/7.801635488463321*234.18157380078924*signal_equiv_peak              # These calibration values lead to a maximum backscatter of the nadir echo of 1.
                                                                            # Therefore, the maximum backsactter can be freely selected via signal_equiv_peak.

u=u1+u2                                                                     # Superposition of scene and nadir echo.


# Plot the data if desired
if (plots):
    plt.figure()
    plt.imshow(20*np.log10(np.abs(u)/float(np.abs(u).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Raw_data_u')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Raw_data_u_befor_dual_focusing')
    plt.show()

    plt.figure()
    plt.imshow(20*np.log10(np.abs(u1)/float(np.abs(u1).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Raw data u1')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Raw_data_u1_befor_dual_focusing')
    plt.show()

    plt.figure()
    plt.imshow(20*np.log10(np.abs(u2)/float(np.abs(u2).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Raw data u2')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Raw_data_u2_befor_dual_focusing')
    plt.show()

print( 'Raw data')
data_values.write('Raw data\n\n')

# Calculate values if desired
if (data):
    data_values.write( 'Peak u1 = '+ str(20.*np.log10(abs(u1).max()))+ ' dB\n')
    data_values.write('Peak u2 = ' + str(20. * np.log10(abs(u2).max())) + ' dB\n')
    data_values.write('Peak u= ' + str(20. * np.log10(abs(u).max())) + ' dB\n\n')

    data_values.write( 'Energy u1= '+ str(10.*np.log10((abs(u1)**2.).sum()))+ ' dB\n')
    data_values.write('Energy u2= ' + str(10. * np.log10((abs(u2) ** 2.).sum())) + ' dB\n')
    data_values.write('Energy u= ' + str(10. * np.log10((abs(u) ** 2.).sum())) + ' dB\n\n\n')

ref_level_raw_data=np.max(abs(u))

# *****
# Filtering adapted to the nadir echo
# *****

# RANGE COMPRESSION

if (hamming_rg):

  n_rg_eff = 2. * np.round((n_rg * B / f_s)/2.)

  h_hamming_rg_F = np.roll(np.concatenate((np.repeat(0., (n_rg - n_rg_eff) / 2.), general_hamming(int(n_rg_eff),hamming_rg_alpha), np.repeat(0., (n_rg - n_rg_eff) / 2.)),axis=0), int(n_rg / 2))
# Hamming window in range


h_rg_desired = np.roll(np.less_equal(np.abs(t_rg-2*R0/c0), tau/2) *np.exp(- j* np.pi * B / tau * (t_rg - 2. * R0 / c0) ** 2.), int(n_rg / 2.))

# impulse response of the range compression filter

temp_var=np.fft.fft(h_rg_desired)*np.conj(np.fft.fft(h_rg_desired))

for k in range(0,n_az):

    t_i_local=t_i[(k+1) % len(t_i)]

    h_rg = np.roll(np.less_equal(np.abs(t_rg - 2 * R0 / c0), tau / 2) * np.exp(j * np.pi * B / tau * (t_rg - 2. * R0 / c0 - t_i_local - tau * np.floor((t_rg - 2. * R0 / c0 + tau / 2 - t_i_local) / tau)) ** 2.), int(n_rg / 2.))

    h_rg_F_ideal = temp_var / np.fft.fft(h_rg)

    if (hamming_rg):
        h_rg_F_ideal = h_rg_F_ideal * h_hamming_rg_F

    u[:, k] = np.fft.ifft(np.fft.fft(u[:, k]) * h_rg_F_ideal)  # range compression (performed in the frequency domain)
    u1[:, k] = np.fft.ifft(np.fft.fft(u1[:, k]) * h_rg_F_ideal)  # range compression (performed in the frequency domain)
    u2[:, k] = np.fft.ifft(np.fft.fft(u2[:, k]) * h_rg_F_ideal)  # range compression (performed in the frequency domain)


# Plot the data if desired
if (plots):
    plt.figure()
    plt.imshow(20*np.log10(np.abs(u)/float(np.abs(u).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Range compressed u (Filtering to nadir echo)')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Range_compressed_u_to_nadir_echo')
    plt.show()

    plt.figure()
    plt.imshow(20*np.log10(np.abs(u1)/float(np.abs(u1).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Range compressed u1 (Filtering to nadir echo)')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Range_compressed_u1_to_nadir_echo')
    plt.show()

    plt.figure()
    plt.imshow(20*np.log10(np.abs(u2)/float(np.abs(u2).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Range compressed u2 (Filtering to nadir echo)')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Range_compressed_u2_to_nadir_echo')
    plt.show()


print( 'Range-compressed data')
data_values.write('Range compressed data to the nadir echo\n\n')

# np.save(folder_name+'/u_dual_rangecompressed_nadir', u)
# np.save(folder_name+'/u1_dual_rangecompressed_nadir', u1)
# np.save(folder_name+'/u2_dual_rangecompressed_nadir', u2)

# Calculate values if desired
if (data):
    data_values.write( 'Peak u1 = '+ str(20.*np.log10(abs(u1).max()))+ ' dB\n')
    data_values.write('Peak u2 = ' + str(20. * np.log10(abs(u2).max())) + ' dB\n')
    data_values.write('Peak u= ' + str(20. * np.log10(abs(u).max())) + ' dB\n\n')

    data_values.write( 'Energy u1= '+ str(10.*np.log10((abs(u1)**2.).sum()))+ ' dB\n')
    data_values.write('Energy u2= ' + str(10. * np.log10((abs(u2) ** 2.).sum())) + ' dB\n')
    data_values.write('Energy u= ' + str(10. * np.log10((abs(u) ** 2.).sum())) + ' dB\n\n\n')

# AZIMUTH COMPRESSION

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
    u[k, :] = np.fft.ifft(np.fft.fft(u[k, :]) * h_az_F)
    u1[k, :] = np.fft.ifft(np.fft.fft(u1[k, :]) * h_az_F)
    u2[k, :] = np.fft.ifft(np.fft.fft(u2[k, :]) * h_az_F)

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



# Display of focused data matched to the nadir (nadir echo focused, signal unfocused/smeared)

# Plot the data if desired
if (plots):
    plt.figure()
    plt.imshow(20*np.log10(np.abs(u)/float(np.abs(u).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Focused data u (Filtering to nadir echo)')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Focused_data_u_to_nadir_echo')
    plt.show()

    plt.figure()
    plt.imshow(20*np.log10(np.abs(u1)/float(np.abs(u1).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Focused data u1 (Filtering to nadir echo)')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Focused_data_u1_to_nadir_echo')
    plt.show()

    plt.figure()
    plt.imshow(20*np.log10(np.abs(u2)/float(np.abs(u2).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Focused data u2 (Filtering to nadir echo)')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Focused_data_u2_to_nadir_echo')
    plt.show()

# np.save(folder_name+'/u_dual_focused_nadir', u)
# np.save(folder_name+'/u1_dual_focused_nadir', u1)
# np.save(folder_name+'/u2_dual_focused_nadir', u2)

print( 'Focused data')
data_values.write('Focused data to the nadir echo\n\n')

# Calculate values if desired
if (data):
    data_values.write( 'Peak u1 = '+ str(20.*np.log10(abs(u1).max()))+ ' dB\n')
    data_values.write('Peak u2 = ' + str(20. * np.log10(abs(u2).max())) + ' dB\n')
    data_values.write('Peak u= ' + str(20. * np.log10(abs(u).max())) + ' dB\n\n')

    data_values.write( 'Energy u1= '+ str(10.*np.log10((abs(u1)**2.).sum()))+ ' dB\n')
    data_values.write('Energy u2= ' + str(10. * np.log10((abs(u2) ** 2.).sum())) + ' dB\n')
    data_values.write('Energy u= ' + str(10. * np.log10((abs(u) ** 2.).sum())) + ' dB\n\n\n')


ref_level=np.max(abs(u))

# np.save(folder_name+'/u_focused_to_nadir_echo', u)
# np.save(folder_name+'/u_focused_to_nadir_echo', u1)
# np.save(folder_name+'/u_focused_to_nadir_echo', u2)

# *****
# Removal of the nadir echo
# *****

max_u=np.max(np.abs(u))

blanked_lines=8                     # Number of lines to be blanked

data_values.write('Nadir echo removal\n\n')
data_values.write('u[8191+4:8191+4'+str(blanked_lines)+', :] = 0.\n')
data_values.write('u1[8191+4:8191+4'+str(blanked_lines)+', :] = 0.\n')
data_values.write('u2[8191+4:8191+4'+str(blanked_lines)+', :] = 0.\n\n\n')


# np.save(folder_name+'/u_focused_nadir', u)
# np.save(folder_name+'/u1_focused_nadir', u1)
# np.save(folder_name+'/u2_focused_nadir', u2)

u[8191+4:8191+4+blanked_lines, :] = 0.
u1[8191+4:8191+4+blanked_lines, :] = 0.
u2[8191+4:8191+4+blanked_lines, :] = 0.

# blanked_parm=nadir_echo_blanking.nadir_echo_blanking(u, u1, u2, 9000, 6)
#
# data_values.write('Blanked dB Value (amp) ='+str(blanked_parm[1])+"dB \n\n")
# data_values.write('Blanked lines ='+str(blanked_parm[0])+"dB \n\n")


# np.save(folder_name+'/u_dual_focused_nadir_after_removal', u)
# np.save(folder_name+'/u1_dual_focused_nadir_after_removal', u1)
# np.save(folder_name+'/u2_dual_focused_nadir_after_removal', u2)

# Plot the data if desired
if (plots):
    plt.figure()
    plt.imshow(20*np.log10(np.abs(u)/float(np.abs(u).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Focused data u2 (Filtering to nadir echo)')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Focused_data_u_to_nadir_echo_nadir_echo_removed')
    plt.show()

print( 'Focused data after nadir removal')
data_values.write('Focused data to the nadir echo after nadir echo removal\n\n')

# Calculate values if desired
if (data):
    data_values.write( 'Peak u1 = '+ str(20.*np.log10(abs(u1).max()))+ ' dB\n')
    data_values.write('Peak u2 = ' + str(20. * np.log10(abs(u2).max())) + ' dB\n')
    data_values.write('Peak u= ' + str(20. * np.log10(abs(u).max())) + ' dB\n\n')

    data_values.write( 'Energy u1= '+ str(10.*np.log10((abs(u1)**2.).sum()))+ ' dB\n')
    data_values.write('Energy u2= ' + str(10. * np.log10((abs(u2) ** 2.).sum())) + ' dB\n')
    data_values.write('Energy u= ' + str(10. * np.log10((abs(u) ** 2.).sum())) + ' dB\n\n\n')

# *****
# # Inverse SAR filtering
#  *****

h_az = np.exp(1j * 4. * np.pi / lamb * h_nad)
t_p = float(R0 * lamb / L / v_g)  # observation time [s]
h_az[np.where(abs(t_az) >= (t_p / 2))] = 0
h_az = np.roll(h_az, int(n_az / 2))
h_az_F = np.fft.fft(h_az)
for k in range(0, n_rg):
 ############# but you have to make a cycle for all ranges
    u[k, :] = np.fft.ifft(np.fft.fft(u[k, :]) / h_az_F)
    u1[k, :] = np.fft.ifft(np.fft.fft(u1[k, :]) / h_az_F)
    u2[k, :] = np.fft.ifft(np.fft.fft(u2[k, :]) / h_az_F)


# Plot the data if desired
if (plots):
    plt.figure()
    plt.imshow(20*np.log10(np.abs(u)/float(np.abs(u).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Defocused data u (Filtering to nadir echo)')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Defocused_data_u_to_nadir_echo')
    plt.show()

    plt.figure()
    plt.imshow(20*np.log10(np.abs(u1)/float(np.abs(u1).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Defocused data u1 (Filtering to nadir echo)')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Defocused_data_u1_to_nadir_echo')
    plt.show()

    plt.figure()
    plt.imshow(20*np.log10(np.abs(u2)/float(np.abs(u2).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Defocused data u2 (Filtering to nadir echo)')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Defocused_data_u2_to_nadir_echo')
    plt.show()

# np.save(folder_name+'/u_dual_rangecompressed_decompressed_nadir', u)
# np.save(folder_name+'/u1_dual_rangecompressed_decompressed_nadir', u1)
# np.save(folder_name+'/u2_dual_rangecompressed_decompressed_nadir', u2)

print( 'Range-compressed data')
data_values.write('Data after de-compression\n\n')

# Calculate values if desired
if (data):
    data_values.write( 'Peak u1 = '+ str(20.*np.log10(abs(u1).max()))+ ' dB\n')
    data_values.write('Peak u2 = ' + str(20. * np.log10(abs(u2).max())) + ' dB\n')
    data_values.write('Peak u= ' + str(20. * np.log10(abs(u).max())) + ' dB\n\n')

    data_values.write( 'Energy u1= '+ str(10.*np.log10((abs(u1)**2.).sum()))+ ' dB\n')
    data_values.write('Energy u2= ' + str(10. * np.log10((abs(u2) ** 2.).sum())) + ' dB\n')
    data_values.write('Energy u= ' + str(10. * np.log10((abs(u) ** 2.).sum())) + ' dB\n\n\n')


# RANGE "DECOMPRESSION"

for k in range(0,n_az):

    t_i_local=t_i[(k+1) % len(t_i)]

    h_rg = np.roll(np.less_equal(np.abs(t_rg - 2 * R0 / c0), tau / 2) * np.exp(j * np.pi * B / tau * (t_rg - 2. * R0 / c0 - t_i_local - tau * np.floor((t_rg - 2. * R0 / c0 + tau / 2 - t_i_local) / tau)) ** 2.), int(n_rg / 2.))

    h_rg_F_ideal = temp_var / np.fft.fft(h_rg)

    if (hamming_rg):
        h_rg_F_ideal = h_rg_F_ideal * h_hamming_rg_F

    n_rg_eff = 2. * round((n_rg * B / f_s) / 2.)

    h_rg_F_ideal_mask = np.roll(np.concatenate((np.repeat(1., (n_rg - n_rg_eff) / 2.), np.repeat(0., n_rg_eff), np.repeat(1., (n_rg - n_rg_eff) / 2.)), axis=0), int(n_rg / 2))  # azimuth window to select the processed Doppler bandwidth


    u[:, k] = np.fft.ifft(np.fft.fft(u[:, k]) * 1. / (h_rg_F_ideal + h_rg_F_ideal_mask) * (1. - h_rg_F_ideal_mask))  # range compression (performed in the frequency domain)
    u1[:, k] = np.fft.ifft(np.fft.fft(u1[:, k]) * 1. / (h_rg_F_ideal + h_rg_F_ideal_mask) * (1. - h_rg_F_ideal_mask))  # range compression (performed in the frequency domain)
    u2[:, k] = np.fft.ifft(np.fft.fft(u2[:, k]) * 1. / (h_rg_F_ideal + h_rg_F_ideal_mask) * (1. - h_rg_F_ideal_mask))  # range compression (performed in the frequency domain)


# Plot the data if desired
if (plots):
    plt.figure()
    plt.imshow(20*np.log10(np.abs(u)/float(np.abs(u).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Data u after range decompression')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Range_decompressed_u_to_nadir_echo')
    plt.show()

    plt.figure()
    plt.imshow(20*np.log10(np.abs(u1)/float(np.abs(u1).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Data u1 after range decompression')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Range_decompressed_u1_to_nadir_echo')
    plt.show()

    plt.figure()
    plt.imshow(20*np.log10(np.abs(u2)/float(np.abs(u2).max())), aspect="auto", cmap='gnuplot2_r', interpolation= 'catrom')
    plt.clim(-100,0)
    plt.colorbar()
    plt.title('Data u2 after range decompression')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Range_decompressed_u2_to_nadir_echo')
    plt.show()

# np.save(folder_name+'/u_dual_rangecompressed_defocused_nadir', u)
# np.save(folder_name+'/u1_dual_rangecompressed_defocused_nadir', u1)
# np.save(folder_name+'/u2_dual_rangecompressed_defocused_nadir', u2)

print( 'Raw data')
data_values.write('Raw data after defocusing\n\n')

# Calculate values if desired
if (data):
    data_values.write( 'Peak u1 = '+ str(20.*np.log10(abs(u1).max()))+ ' dB\n')
    data_values.write('Peak u2 = ' + str(20. * np.log10(abs(u2).max())) + ' dB\n')
    data_values.write('Peak u= ' + str(20. * np.log10(abs(u).max())) + ' dB\n\n')

    data_values.write( 'Energy u1= '+ str(10.*np.log10((abs(u1)**2.).sum()))+ ' dB\n')
    data_values.write('Energy u2= ' + str(10. * np.log10((abs(u2) ** 2.).sum())) + ' dB\n')
    data_values.write('Energy u= ' + str(10. * np.log10((abs(u) ** 2.).sum())) + ' dB\n\n\n')



# *****
# Filtering adapted to the signal
# *****

# RANGE COMPRESSION

for k in range(0,n_az):

    t_i_local=t_i[k % len(t_i)]

    h_rg = np.roll(np.less_equal(np.abs(t_rg - 2 * R0 / c0), tau / 2) * np.exp(j * np.pi * B / tau * (t_rg - 2. * R0 / c0 - t_i_local - tau * np.floor((t_rg - 2. * R0 / c0 + tau / 2 - t_i_local) / tau)) ** 2.), int(n_rg / 2.))

    h_rg_F_ideal = temp_var / np.fft.fft(h_rg)

    if (hamming_rg):
        h_rg_F_ideal = h_rg_F_ideal * h_hamming_rg_F


    u[:, k] = np.fft.ifft(np.fft.fft(u[:, k]) * h_rg_F_ideal)       # range compression (performed in the frequency domain)
    u1[:, k] = np.fft.ifft(np.fft.fft(u1[:, k]) * h_rg_F_ideal)     # range compression (performed in the frequency domain)
    u2[:, k] = np.fft.ifft(np.fft.fft(u2[:, k]) * h_rg_F_ideal)     # range compression (performed in the frequency domain)


# Plot the data if desired
if (plots):
    plt.figure()
    plt.imshow(20 * np.log10(np.abs(u) / float(np.abs(u).max())), aspect="auto", cmap='gnuplot2_r', interpolation='catrom')
    plt.clim(-100, 0)
    plt.colorbar()
    plt.title('u range compressed to the signal')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Range_compressed data_u_to_the_signal_after_nadir_echo_removal')
    plt.show()

    plt.figure()
    plt.imshow(20 * np.log10(np.abs(u1) / float(np.abs(u1).max())), aspect="auto", cmap='gnuplot2_r',
               interpolation='catrom')
    plt.clim(-100, 0)
    plt.colorbar()
    plt.title('u1 range compressed to the signal')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Range_compressed data_u1_to_the_signal_after_nadir_echo_removal')
    plt.show()

    plt.figure()
    plt.imshow(20 * np.log10(np.abs(u2) / float(np.abs(u2).max())), aspect="auto", cmap='gnuplot2_r',
               interpolation='catrom')
    plt.clim(-100, 0)
    plt.colorbar()
    plt.title('u2 range compressed to the signal')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Range_compressed data_u2_to_the_signal_after_nadir_echo_removal')
    plt.show()

# np.save(folder_name+'/u_dual_rangecompressed_signal', u)
# np.save(folder_name+'/u1_dual_rangecompressed_signal', u1)
# np.save(folder_name+'/u2_dual_rangecompressed_signal', u2)

print( 'Range-compressed data')
data_values.write('Range compressed data to the signal after nadir echo removal \n\n')

# Calculate values if desired
if (data):
    data_values.write( 'Peak u1 = '+ str(20.*np.log10(abs(u1).max()))+ ' dB\n')
    data_values.write('Peak u2 = ' + str(20. * np.log10(abs(u2).max())) + ' dB\n')
    data_values.write('Peak u= ' + str(20. * np.log10(abs(u).max())) + ' dB\n\n')

    data_values.write( 'Energy u1= '+ str(10.*np.log10((abs(u1)**2.).sum()))+ ' dB\n')
    data_values.write('Energy u2= ' + str(10. * np.log10((abs(u2) ** 2.).sum())) + ' dB\n')
    data_values.write('Energy u= ' + str(10. * np.log10((abs(u) ** 2.).sum())) + ' dB\n\n\n')

# AZIMUTH COMPRESSION

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
    u[k, :] = np.fft.fft(u[k, :]) * h_az_F
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
    u[k, :] = np.fft.ifft(u[k, :] * az_window_F)
    u1[k, :] = np.fft.ifft(u1[k, :] * az_window_F)
    u2[k, :] = np.fft.ifft(u2[k, :] * az_window_F)



# Plot the data if desired
if (plots):
    plt.figure()
    plt.imshow(20 * np.log10(np.abs(u) / float(np.abs(u).max())), aspect="auto", cmap='gnuplot2_r',
               interpolation='catrom')
    plt.clim(-100, 0)
    plt.colorbar()
    plt.title('u focused to the signal')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Focused data_u_to_the_signal_after_nadir_echo_removal')
    plt.show()

    plt.figure()
    plt.imshow(20 * np.log10(np.abs(u1) / float(np.abs(u1).max())), aspect="auto", cmap='gnuplot2_r',
               interpolation='catrom')
    plt.clim(-100, 0)
    plt.colorbar()
    plt.title('u1 focused to the signal')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Focused data_u1_to_the_signal_after_nadir_echo_removal')
    plt.show()

    plt.figure()
    plt.imshow(20 * np.log10(np.abs(u2) / float(np.abs(u2).max())), aspect="auto", cmap='gnuplot2_r',
               interpolation='catrom')
    plt.clim(-100, 0)
    plt.colorbar()
    plt.title('u2 focused to the signal')
    plt.ylabel('Range samples')
    plt.xlabel('Azimuth samples')
    plt.savefig(folder_name+'\Focused data_u2_to_the_signal_after_nadir_echo_removal')
    plt.show()

print( 'Focused data')
data_values.write('Focused data to the signal after nadir echo removal \n\n')

# Calculate values if desired
if (data):
    data_values.write( 'Peak u1 = '+ str(20.*np.log10(abs(u1).max()))+ ' dB\n')
    data_values.write('Peak u2 = ' + str(20. * np.log10(abs(u2).max())) + ' dB\n')
    data_values.write('Peak u= ' + str(20. * np.log10(abs(u).max())) + ' dB\n\n')

    data_values.write( 'Energy u1= '+ str(10.*np.log10((abs(u1)**2.).sum()))+ ' dB\n')
    data_values.write('Energy u2= ' + str(10. * np.log10((abs(u2) ** 2.).sum())) + ' dB\n')
    data_values.write('Energy u= ' + str(10. * np.log10((abs(u) ** 2.).sum())) + ' dB\n\n\n')


data_values.close()

# Save processed data after dual focusing u, u1, u2
np.save(folder_name+'/u_dual_focused', u)
np.save(folder_name+'/u1_dual_focused', u1)
np.save(folder_name+'/u2_dual_focused', u2)