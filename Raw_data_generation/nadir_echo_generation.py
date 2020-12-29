'''
##########
# author Fabian Glatz
# Description:
# With this code, a signal (complex data matrix) of a realistic nadir echo can be generated.
# The model is based on the parametric nadir echo model from Michelangelo and Maxwell
##########
'''

# Imports
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal.windows import general_hamming
import Parameters
import os
from datetime import datetime

# Paramters

param = Parameters.Parameters()         # Load dia parameters from parameters class

n_rg = param.n_rg
n_az = param.n_az

t_rg = param.t_rg
t_az = param.t_az

PRF = param.PRF


t_max = t_rg[8196]              # Sample in range direction where the nadir echo has the maximum backscatter
tau_rise = 28.03*10**-9         # rise time of the nadir echo
tau_fall = 60.97*10**-9         # fall time of the nadir echo (max to half of the intensity value)
alpha = 2.4                     # parameter to adjust the shape of the rising edge
beta = 0.33                     # parameter to adjust the shape of the falling edge
sigma_max = 1                   # value of the maximum backscatter


def fun (t_max,tau_rise,tau_fall,alpha,beta,sigma_max, t_rg):
    if t_rg < t_max-tau_rise:
        return 0
    elif t_max-tau_rise <= t_rg < t_max-tau_rise/2:
        return sigma_max*(2**(alpha-1))*((t_rg-t_max)/tau_rise+1)**alpha
    elif t_max-tau_rise/2 <= t_rg < t_max:
        return sigma_max*(1-(2**(alpha-1))*((t_max-t_rg)/tau_rise)**alpha)
    elif t_max <= t_rg < t_max + tau_fall:
        return sigma_max * ((1 - beta) * ((t_rg - t_max) / tau_fall) ** 3 + (beta - 3 / 2) * ((t_rg - t_max) / tau_fall) ** 2 + 1)
    elif t_max + tau_fall <= t_rg:
        return sigma_max/(2*((t_rg-t_max)/tau_fall)**(2*beta))


nadir_echo = np.zeros((len(t_rg)))
for k in range(len(t_rg)):
    nadir_echo[k]=fun(t_max, tau_rise, tau_fall, alpha, beta, sigma_max, t_rg[k])

# Other way to define a piecewise function
#nadir_echo=np.piecewise(t_rg, [t_rg < t_max-tau_rise,  t_max-tau_rise <= t_rg, t_max-tau_rise/2 <= t_rg, t_max <= t_rg, t_max+tau_fall <= t_rg], [lambda t_rg: 0, lambda t_rg: sigma_max*(2**(alpha-1))*((t_rg-t_max)/tau_rise+1)**alpha, lambda t_rg: sigma_max*(1-(2**(alpha-1))*((t_max-t_rg)/tau_rise)**alpha), lambda t_rg: sigma_max*((1-beta)*((t_rg-t_max)/tau_fall)**3+(beta-3/2)*((t_rg-t_max)/tau_fall)**2+1) , lambda t_rg: sigma_max/(2*((t_rg-t_max)/tau_fall)**(2*beta))])



# Plot cutplot of the nadir echo in range direction

l=len(t_rg)

plt.plot(t_rg, nadir_echo)
plt.title('Cutplot of the nadir echo in range direction')
plt.xlabel('range time [s]')
plt.ylabel('Amplitude')
plt.show()

plt.plot(t_rg[int(l/2-10):int(l/2+100)],nadir_echo[int(l/2-10):int(l/2+100)])
plt.title('Enlarged viewing cutplot of the nadir echo in range direction')
plt.xlabel('range time [s]')
plt.ylabel('Amplitude')
plt.show()


# Extend the nadir echo shape to all azimuth samples (The tilt due to the curvature of the earth is neglected)

nadir_echo_1 = np.ones((n_rg, n_az))                    # Generate a matrix with ones n_rg x n_az

nadir_echo=np.array([nadir_echo])                       # change the line vector to a vector with 1 x n_rg just to make it possible to transpose the vector

backscatter_nadir_echo = nadir_echo_1 * nadir_echo.T    # do a vector multiplication to extend the shape of the nadir echo to all azimuth samples


# Plot the 2D intensity backscatter
plt.imshow(np.abs(backscatter_nadir_echo),cmap="gray",vmin=0,vmax=2*np.mean(np.abs(backscatter_nadir_echo)))
plt.title('Backscatter of the nadir echo')
plt.xlabel('azimuth samples')
plt.ylabel('range samples')
plt.show()

# Generate speckle matrix to add the speckle effect to the nadir echo backscatter

speckle =  1./np.sqrt(2)*(np.random.randn(n_rg, n_az)+1j*np.random.randn(n_rg, n_az)) # generate a complex speckle matrix

nadir_echo_1 = backscatter_nadir_echo*speckle   # muliply the speckle vector with the backscatter of the nadir echo


# Plot the 2D intensity backscatter after adding speckle
plt.imshow(np.abs(nadir_echo_1),cmap="gray",vmin=0,vmax=2*np.mean(np.abs(nadir_echo_1)))
plt.title('Backscatter of the nadir echo with speckle')
plt.xlabel('azimuth samples')
plt.ylabel('range samples')
plt.show()

# Plot enlarged shape of the nadir echos after adding speckle
plt.plot(t_rg[int(l/2-10):int(l/2+100)], np.abs(nadir_echo_1[int(l/2-10):int(l/2+100),1000]))
plt.title('Enlarged viewing: Cutplot of the nadir echo in range \n direction after adding speckle')
plt.xlabel('range time [s]')
plt.ylabel('Intensity')
plt.show()

# Autocorrelation function of the nadir echo

R_nad = np.exp(-(t_az/(0.46*10**-3))**2)  # Formula from Michelangelo and Maxwell

# Plot autocorrelation function of the nadir echo
plt.plot(t_az, np.real(R_nad))
plt.title('Autocorrelation function (not oversampled)')
plt.xlabel('azimuth time [s]')
plt.ylabel('Amplitude')
plt.show()

plt.plot(t_az[int(len(R_nad)/2)-10:int(len(R_nad)/2)+11],np.real(R_nad[int(len(R_nad)/2)-10:int(len(R_nad)/2)+11]))
plt.title('Enlarged viewing: Autocorrelation function (not oversampled)')
plt.xlabel('azimuth time [s]')
plt.ylabel('Amplitude')
plt.show()

############### Oversampling of the autocorrelation function ############
# Just for better imaging
#f_az = np.roll((np.arange(n_az) - n_az / 2.) * PRF / n_az, int(n_az / 2.))

# R_nad_F=np.fft.fft(R_nad) # fft of the autocorrelation function
#
# R_nad_F_zero=np.concatenate([R_nad_F[0:n_az//2], np.zeros(31*n_az,'complex'), R_nad_F[n_az//2:n_az]])
#
# R_nad_over=np.fft.ifft(R_nad_F_zero)*32
#
# t_az_zero = (np.arange(len(R_nad_over)) - len(R_nad_over) / 2) / (PRF*32)
#
# plt.plot(t_az_zero, np.real(R_nad_over))
# plt.title('Autocorrelation function (oversampled)')
# plt.xlabel('azimuth time [s]')
# plt.ylabel('Amplitude')
# plt.show()
#
# plt.plot(t_az_zero[int(len(R_nad_over)/2)-310:int(len(R_nad_over)/2)+311], np.real(R_nad_over[int(len(R_nad_over)/2)-310:int(len(R_nad_over)/2)+311]))
# plt.title('Enlarged viewing: Autocorrelation function (oversampled)')
# plt.xlabel('azimuth time [s]')
# plt.ylabel('Amplitude')
# plt.show()

#########################################################################

# Calculate the power spectral density of nadir echos autocorrelation function

P_desired = np.fft.fft(np.fft.fftshift(R_nad))

## Plot the power spectral density
# f_az = (np.arange(n_az) - n_az / 2.) * PRF / n_az
#
# plt.plot(f_az, np.real(np.fft.fftshift(P_desired)))
# plt.title('P_desired')
# plt.xlabel('frequency [Hz]')
# plt.ylabel('Amplitude')
# plt.show()

nadir_echo_2 = np.fft.ifft(np.fft.fft(nadir_echo_1, axis = 1)*np.sqrt(P_desired), axis =1)  # Add the autocorrelation function in the frequency domain.


# Sqrt to get the amplitude values
nadir_echo_2=np.sqrt(nadir_echo_2)          # Transfer the intensity values to amplitude values

# Save the generated nadir echo
np.save('nadir_echo', nadir_echo_2)


############### Oversampling of the autocorrelation function ############

# line_az_F=np.fft.fft(nadir_echo_2[10000,:])
#
# line_az_F_zero=np.concatenate([line_az_F[0:n_az//2], np.zeros(31*n_az,'complex'), line_az_F[n_az//2:n_az]])
#
# line_az_over=np.fft.ifft(line_az_F_zero)*32
#
# line_az_over_auto = np.correlate(line_az_over, line_az_over, mode='same')
#
# plt.plot(t_az_zero[int(len(line_az_over_auto)/2)-310:int(len(line_az_over_auto)/2)+311], np.real(line_az_over_auto[int(len(line_az_over_auto)/2)-310:int(len(line_az_over_auto)/2)+311]))
# plt.title('Autocorrelation of one single line in azimuth direction \n after adding the typical autocorrelation of the nadir echo')
# plt.xlabel('azimuth time [s]')
# plt.ylabel('Amplitude')
# plt.show()
#
# line_az_F=np.fft.fft(nadir_echo_1[10000,:])
#
# line_az_F_zero=np.concatenate([line_az_F[0:n_az//2], np.zeros(31*n_az,'complex'), line_az_F[n_az//2:n_az]])
#
# line_az_over=np.fft.ifft(line_az_F_zero)*32
#
# line_az_over_auto = np.correlate(line_az_over, line_az_over, mode='same')
#
# plt.plot(t_az_zero[int(len(line_az_over_auto)/2)-310:int(len(line_az_over_auto)/2)+311], np.real(line_az_over_auto[int(len(line_az_over_auto)/2)-310:int(len(line_az_over_auto)/2)+311]))
# plt.title('Autocorrelation of one single line in azimuth direction \n without adding the autocorrelation of the nadir echo')
# plt.xlabel('azimuth time [s]')
# plt.ylabel('Amplitude')
# plt.show()
#



# Estimation of the autocorrelation function

# r = np.zeros(n_az, dtype=complex)
#
# for k in range(0, 20):
#     r[k] = np.mean(np.conj(nadir_echo_2[10000,:]) * np.roll(nadir_echo_2[10000,:], -k)) / np.mean(abs(nadir_echo_2[10000,:]) ** 2.)
#
# for k in range(1, 21):
#     r[n_az-k] = np.mean(np.conj(nadir_echo_2[10000,:]) * np.roll(nadir_echo_2[10000,:], k)) / np.mean(abs(nadir_echo_2[10000,:]) ** 2.)
#
# test=np.fft.fftshift(r)
#
# plt.plot(test[4096-10:4096+11])
# plt.show()

#for k=0,20 do r(k)=mean(conj(u_dem)*shift(u_dem,0,-k))/mean(abs(u_dem)^2.) IDL

#for k=1,21 do r(n-k)=mean(conj(u_dem)*shift(u_dem,0,k))/mean(abs(u_dem)^2.) IDL