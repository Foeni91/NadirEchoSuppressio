"""
# author Fabian Glatz
##################################################################
# Description:
# This code shows the generation of different probability density functions.
##################################################################
"""
import numpy as np
import matplotlib.pyplot as plt
import Parameters

param = Parameters.Parameters()


path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_Set_2_normal/Up_down_Chirp_no_RCM_CR_5050_n_az_2_13_20200904_174912" #Put here the folder name wher the data is located and where the images are store

u=np.load(path_temp+'/u2_single_focused.npy')

path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_Set_2_L_Band/Up_down_Chirp_no_RCM_CR_5050_n_az_2_13_20200913_141356"
u1=np.load(path_temp+'/u2_single_focused.npy')

path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_Set_2_normal/Shifted_Chirp_no_RCM_CR_5050_n_az_2_13_20200912_184619"
u2=np.load(path_temp+'/u2_single_focused.npy')

path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_Set_2_L_Band/Shifted_Chirp_no_RCM_CR_5050_n_az_2_13_20200913_125457"
u3=np.load(path_temp+'/u2_single_focused.npy')

u_single=u[6144:10240,2048:6144]
u1_single=u1[6144:10240,2048:6144]
u2_single=u2[6144:10240,2048:6144]
u3_single=u3[6144:10240,2048:6144]

# region Probability density function of the log-intensity


# Compute the histograms
bin_num=200

hist_int_log_u_single, bin_edge_log_u_single = np.histogram(10*np.log10((np.abs(u_single) ** 2)), bins=bin_num)
hist_int_log_u1_single, bin_edge_log_u1_single = np.histogram(10*np.log10((np.abs(u1_single) ** 2)), bins=bin_num)
hist_int_log_u2_single, bin_edge_log_u2_single = np.histogram(10*np.log10((np.abs(u2_single) ** 2)), bins=bin_num)
hist_int_log_u3_single, bin_edge_log_u3_single = np.histogram(10*np.log10((np.abs(u3_single) ** 2)), bins=bin_num)

# Display the probability density function

# City
plt.figure()
plt.semilogy(bin_edge_log_u_single[:-1], hist_int_log_u_single/(u_single.size)/(bin_edge_log_u_single[2]-bin_edge_log_u_single[1]),'b-', label="Up-Down X-Band")
plt.semilogy(bin_edge_log_u1_single[:-1], hist_int_log_u1_single/(u1_single.size)/(bin_edge_log_u1_single[2]-bin_edge_log_u1_single[1]),'g-', label="Up-Down L-Band")
plt.semilogy(bin_edge_log_u2_single[:-1], hist_int_log_u2_single/(u2_single.size)/(bin_edge_log_u2_single[2]-bin_edge_log_u2_single[1]),'r-', label="Shifted X-Band")
plt.semilogy(bin_edge_log_u3_single[:-1], hist_int_log_u3_single/(u3_single.size)/(bin_edge_log_u3_single[2]-bin_edge_log_u3_single[1]),'-',color='darkviolet', label="Shifted L-Band")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)




u=np.load(path_temp+'/u_dual_focused.npy')
u1=np.load(path_temp+'/u1_dual_focused.npy')
u2=np.load(path_temp+'/u2_dual_focused.npy')

u_dual=u[6144+8192:10240+8192,2048:6144]
u1_dual=u1[6144+8192:10240+8192,2048:6144]
u2_dual=u2[6144+8192:10240+8192,2048:6144]

# region Probability density function of the log-intensity


# Compute the histograms
bin_num=200

hist_int_log_u_dual, bin_edge_log_u_dual = np.histogram(10*np.log10((np.abs(u_dual) ** 2)), bins=bin_num)
hist_int_log_u1_dual, bin_edge_log_u1_dual = np.histogram(10*np.log10((np.abs(u1_dual) ** 2)), bins=bin_num)
hist_int_log_u2_dual, bin_edge_log_u2_dual = np.histogram(10*np.log10((np.abs(u2_dual) ** 2)), bins=bin_num)


# Display the probability density function

# City
plt.figure()
plt.semilogy(bin_edge_log_u_dual[:-1], hist_int_log_u_dual/(u_dual.size)/(bin_edge_log_u_dual[2]-bin_edge_log_u_dual[1]),'b-', label="u_dual")
plt.semilogy(bin_edge_log_u1_dual[:-1], hist_int_log_u1_dual/(u1_dual.size)/(bin_edge_log_u1_dual[2]-bin_edge_log_u1_dual[1]),'g-', label="u1_dual")
plt.semilogy(bin_edge_log_u2_dual[:-1], hist_int_log_u2_dual/(u2_dual.size)/(bin_edge_log_u2_dual[2]-bin_edge_log_u2_dual[1]),'r-', label="u2_dual")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)
plt.savefig(path_temp+"/u_u1_u2_dual_pdf_log")


# Long short

u1_long=np.load('Up_down_Chirp_no_RCM_CR_5050_n_az_2_13_20200904_174912/u1_single_focused.npy')
u1_short=np.load('Up_down_Chirp_no_RCM_CR_2500_n_az_2_13_20200905_142132/u1_single_focused.npy')
u1_single_long=u1_long[6144:10240,2048:6144]
u1_single_short=u1_short[6144:10240,2048:6144]

bin_num=200
hist_int_log_u1_single_short, bin_edge_log_u1_single_short = np.histogram(10*np.log10((np.abs(u1_single_short) ** 2)), bins=bin_num)
hist_int_log_u1_single_long, bin_edge_log_u1_single_long = np.histogram(10*np.log10((np.abs(u1_single_long) ** 2)), bins=bin_num)

plt.figure()
plt.semilogy(bin_edge_log_u1_single_short[:-1], hist_int_log_u1_single_short/(u1_single_short.size)/(bin_edge_log_u1_single_short[2]-bin_edge_log_u1_single_short[1]),'b-', label="u_single_short")
plt.semilogy(bin_edge_log_u1_single_long[:-1], hist_int_log_u1_single_long/(u1_single_long.size)/(bin_edge_log_u1_single_long[2]-bin_edge_log_u1_single_long[1]),'r-', label="u_single_long")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)


bin_num=200
hist_int_log_u2_single_short, bin_edge_log_u2_single_short = np.histogram(10*np.log10((np.abs(u2_single_short) ** 2)), bins=bin_num)
hist_int_log_u2_single_long, bin_edge_log_u2_single_long = np.histogram(10*np.log10((np.abs(u2_single_long) ** 2)), bins=bin_num)

plt.figure()
plt.semilogy(bin_edge_log_u2_single_short[:-1], hist_int_log_u2_single_short/(u2_single_short.size)/(bin_edge_log_u2_single_short[2]-bin_edge_log_u2_single_short[1]),'b-', label="u2_single_short")
plt.semilogy(bin_edge_log_u2_single_long[:-1], hist_int_log_u2_single_long/(u2_single_long.size)/(bin_edge_log_u2_single_long[2]-bin_edge_log_u2_single_long[1]),'r-', label="u2_single_long")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)



#Nadir echo noise level#

# Conv SAR

path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_Set_2_normal/Conv_no_RCM_CR_5050_n_az_2_13_20200918_172806" #Put here the folder name wher the data is located and where the images are store

u=np.load(path_temp+'/u2_single_focused.npy')

u_single=u[8190:10240,2048:6144]


# Compute the histograms
bin_num=200

hist_int_log_u_single, bin_edge_log_u_single = np.histogram(10*np.log10((np.abs(u_single) ** 2)), bins=bin_num)

plt.figure()
plt.semilogy(bin_edge_log_u_single[:-1], hist_int_log_u_single/(u_single.size)/(bin_edge_log_u_single[2]-bin_edge_log_u_single[1]),'g-', label="Conv. Nadir Echo")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.xlim(-90,15)
plt.ylim(10**-7,0)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)


# NESZ

u1=np.load(path_temp+'/u1_single_focused.npy')
u1_single=u1[7800:8500,2500:4200]

hist_int_log_u1_single, bin_edge_log_u1_single = np.histogram(10*np.log10((np.abs(u1_single) ** 2)), bins=bin_num)

plt.figure()
plt.semilogy(bin_edge_log_u1_single[:-1], hist_int_log_u1_single/(u1_single.size)/(bin_edge_log_u1_single[2]-bin_edge_log_u1_single[1]),'r-', label="NESZ")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)


#Smeard nadir echo

path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_Set_2_normal/Up_down_Chirp_no_RCM_CR_5050_n_az_2_13_20200904_174912" #Put here the folder name wher the data is located and where the images are store

u=np.load(path_temp+'/u2_single_focused.npy')

u_single=u[6144:10240,2048:6144]


# Compute the histograms
bin_num=200

hist_int_log_u_single, bin_edge_log_u_single = np.histogram(10*np.log10((np.abs(u_single) ** 2)), bins=bin_num)

plt.figure()
plt.semilogy(bin_edge_log_u_single[:-1], hist_int_log_u_single/(u_single.size)/(bin_edge_log_u_single[2]-bin_edge_log_u_single[1]),'g-', label="NENE")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)



############### Diff beta with and without blanking ###################


path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/Data_table/5050/Up_down/Up_down_Chirp_no_RCM_CR_5050_n_az_2_13_20201101_142414" #Put here the folder name wher the data is located and where the images are store

u=np.load(path_temp+'/u2_single_focused.npy')
u1=np.load(path_temp+'/u2_dual_focused.npy')

path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/Data_table/5050/Up_down/Up_down_Chirp_no_RCM_CR_5050_n_az_2_13_20201101_153306_beta_07_1"
u2=np.load(path_temp+'/u2_single_focused.npy')
u3=np.load(path_temp+'/u2_dual_focused.npy')

u_single=u[6144:10240,2048:6144]
u1_single=u1[6144:10240,2048:6144]
u2_single=u2[6144:10240,2048:6144]
u3_single=u3[6144:10240,2048:6144]

# region Probability density function of the log-intensity


# Compute the histograms
bin_num=200

hist_int_log_u_single, bin_edge_log_u_single = np.histogram(10*np.log10((np.abs(u_single) ** 2)), bins=bin_num)
hist_int_log_u1_single, bin_edge_log_u1_single = np.histogram(10*np.log10((np.abs(u1_single) ** 2)), bins=bin_num)
hist_int_log_u2_single, bin_edge_log_u2_single = np.histogram(10*np.log10((np.abs(u2_single) ** 2)), bins=bin_num)
hist_int_log_u3_single, bin_edge_log_u3_single = np.histogram(10*np.log10((np.abs(u3_single) ** 2)), bins=bin_num)

# Display the probability density function

# City
plt.figure()
plt.semilogy(bin_edge_log_u_single[:-1], hist_int_log_u_single/(u_single.size)/(bin_edge_log_u_single[2]-bin_edge_log_u_single[1]),'r-', label=r"$\beta = 0.33$")
plt.semilogy(bin_edge_log_u1_single[:-1], hist_int_log_u1_single/(u1_single.size)/(bin_edge_log_u1_single[2]-bin_edge_log_u1_single[1]),'-',color='maroon', label=r"$\beta = 0.33 (dual focused)$")
plt.semilogy(bin_edge_log_u2_single[:-1], hist_int_log_u2_single/(u2_single.size)/(bin_edge_log_u2_single[2]-bin_edge_log_u2_single[1]),'g-', label=r"$\beta = 0.7$")
plt.semilogy(bin_edge_log_u3_single[:-1], hist_int_log_u3_single/(u3_single.size)/(bin_edge_log_u3_single[2]-bin_edge_log_u3_single[1]),'-',color='purple', label=r"$\beta = 0.7 (dual focused)$")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)




######### Up_down and shifted (not propper smeared) ####

path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_Set_2_normal/Same_calibration/Up_down_Chirp_no_RCM_CR_5050_n_az_2_13_20201020_151336" #Put here the folder name wher the data is located and where the images are store

u=np.load(path_temp+'/u2_single_focused.npy')
#u1=np.load(path_temp+'/u2_dual_focused.npy')

path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_Set_2_normal/Same_calibration/Shifted_Chirp_no_RCM_CR_5050_n_az_2_13_20201020_155021"
u2=np.load(path_temp+'/u2_single_focused.npy')
#u3=np.load(path_temp+'/u2_dual_focused.npy')

u_single=u[6144:10240,2048:6144]
#u1_single=u1[6144:10240,2048:6144]
u2_single=u2[6144:10240,2048:6144]
#u3_single=u3[6144:10240,2048:6144]

# region Probability density function of the log-intensity


# Compute the histograms
bin_num=200

hist_int_log_u_single, bin_edge_log_u_single = np.histogram(10*np.log10((np.abs(u_single) ** 2)), bins=bin_num)
#hist_int_log_u1_single, bin_edge_log_u1_single = np.histogram(10*np.log10((np.abs(u1_single) ** 2)), bins=bin_num)
hist_int_log_u2_single, bin_edge_log_u2_single = np.histogram(10*np.log10((np.abs(u2_single) ** 2)), bins=bin_num)
#hist_int_log_u3_single, bin_edge_log_u3_single = np.histogram(10*np.log10((np.abs(u3_single) ** 2)), bins=bin_num)

# Display the probability density function

# City
plt.figure()
plt.semilogy(bin_edge_log_u_single[:-1], hist_int_log_u_single/(u_single.size)/(bin_edge_log_u_single[2]-bin_edge_log_u_single[1]),'b-', label=r"Up- and down-chirps")
#plt.semilogy(bin_edge_log_u1_single[:-1], hist_int_log_u1_single/(u1_single.size)/(bin_edge_log_u1_single[2]-bin_edge_log_u1_single[1]),'-',color='maroon', label=r"$\beta = 0.33 (dual focused)$")
plt.semilogy(bin_edge_log_u2_single[:-1], hist_int_log_u2_single/(u2_single.size)/(bin_edge_log_u2_single[2]-bin_edge_log_u2_single[1]),'-',color='purple', label=r"Cyclically-shifted chirps")
#plt.semilogy(bin_edge_log_u3_single[:-1], hist_int_log_u3_single/(u3_single.size)/(bin_edge_log_u3_single[2]-bin_edge_log_u3_single[1]),'-',color='lime', label=r"$\beta = 0.7 (dual focused)$")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)


######### Up_down and shifted propper smeard K=10 ####

path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_Set_2_normal/Up_down_Chirp_no_RCM_CR_5050_n_az_2_13_20200904_174912" #Put here the folder name wher the data is located and where the images are store

u=np.load(path_temp+'/u2_single_focused.npy')
#u1=np.load(path_temp+'/u2_dual_focused.npy')

path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_Set_2_normal/Other_K/Shifted_Chirp_no_RCM_CR_5050_n_az_2_13_20200919_093957_K_10"
u2=np.load(path_temp+'/u2_single_focused.npy')
#u3=np.load(path_temp+'/u2_dual_focused.npy')

u_single=u[6144:10240,2048:6144]
#u1_single=u1[6144:10240,2048:6144]
u2_single=u2[6144:10240,2048:6144]
#u3_single=u3[6144:10240,2048:6144]

# region Probability density function of the log-intensity


# Compute the histograms
bin_num=200

hist_int_log_u_single, bin_edge_log_u_single = np.histogram(10*np.log10((np.abs(u_single) ** 2)), bins=bin_num)
#hist_int_log_u1_single, bin_edge_log_u1_single = np.histogram(10*np.log10((np.abs(u1_single) ** 2)), bins=bin_num)
hist_int_log_u2_single, bin_edge_log_u2_single = np.histogram(10*np.log10((np.abs(u2_single) ** 2)), bins=bin_num)
#hist_int_log_u3_single, bin_edge_log_u3_single = np.histogram(10*np.log10((np.abs(u3_single) ** 2)), bins=bin_num)

# Display the probability density function

# City
plt.figure()
plt.semilogy(bin_edge_log_u_single[:-1], hist_int_log_u_single/(u_single.size)/(bin_edge_log_u_single[2]-bin_edge_log_u_single[1]),'b-', label=r"Up-down chirp")
#plt.semilogy(bin_edge_log_u1_single[:-1], hist_int_log_u1_single/(u1_single.size)/(bin_edge_log_u1_single[2]-bin_edge_log_u1_single[1]),'-',color='maroon', label=r"$\beta = 0.33 (dual focused)$")
plt.semilogy(bin_edge_log_u2_single[:-1], hist_int_log_u2_single/(u2_single.size)/(bin_edge_log_u2_single[2]-bin_edge_log_u2_single[1]),'-',color='purple', label=r"Shifted chirp K=10")
#plt.semilogy(bin_edge_log_u3_single[:-1], hist_int_log_u3_single/(u3_single.size)/(bin_edge_log_u3_single[2]-bin_edge_log_u3_single[1]),'-',color='lime', label=r"$\beta = 0.7 (dual focused)$")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)



###### PDF Noise exponential distributed ######


noise=1/np.sqrt(2)*(np.random.randn(1000, 1000)+1j*np.random.randn(1000, 1000))

bin_num=200

u_single=noise

hist_int_log_u_single, bin_edge_log_u_single = np.histogram(10*np.log10((np.abs(u_single) ** 2)), bins=bin_num)
hist_int_city, bin_edge_city = np.histogram((np.abs(u_single) ** 2), bins=bin_num)

plt.figure()
plt.semilogy(bin_edge_log_u_single[:-1], hist_int_log_u_single/(u_single.size)/(bin_edge_log_u_single[2]-bin_edge_log_u_single[1]),'b-', label=r"Noise")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)



# Display the probability density function

# City
plt.figure()
plt.plot(bin_edge_city[:-1], hist_int_city/(noise.size)/(bin_edge_city[2]-bin_edge_city[1]),'b-', label="Noise")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [linear]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)



############### Diff pulse width ###################


path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_f_s_110MHz_tau_25us_n_az_8192_n_rg_16384/Up_down_Chirp_no_RCM_CR_2500_n_az_2_13_20200926_144602" #Put here the folder name wher the data is located and where the images are store
u=np.load(path_temp+'/u2_single_focused.npy')

path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_Set_2_normal/Up_down_Chirp_no_RCM_CR_5050_n_az_2_13_20200830_180657_200_lines_blanked"
u1=np.load(path_temp+'/u2_single_focused.npy')

path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_f_s_110MHz_tau_101us_n_az_8192_n_rg_32768/nadir_echo_long_generated/Up_down_Chirp_no_RCM_CR_10100_n_az_2_13_20201030_101812"
u2=np.load(path_temp+'/u2_single_focused.npy')


u_single=u[6144:10240,2048:6144]
u1_single=u1[6144:10240,2048:6144]
u2_single=u2[6144+8192:10240+8192,2048:6144]


# region Probability density function of the log-intensity


# Compute the histograms
bin_num=200

hist_int_log_u_single, bin_edge_log_u_single = np.histogram(10*np.log10((np.abs(u_single) ** 2)), bins=bin_num)
hist_int_log_u1_single, bin_edge_log_u1_single = np.histogram(10*np.log10((np.abs(u1_single) ** 2)), bins=bin_num)
hist_int_log_u2_single, bin_edge_log_u2_single = np.histogram(10*np.log10((np.abs(u2_single) ** 2)), bins=bin_num)


# Display the probability density function

# City
plt.figure()
plt.semilogy(bin_edge_log_u_single[:-1], hist_int_log_u_single/(u_single.size)/(bin_edge_log_u_single[2]-bin_edge_log_u_single[1]),'b-', label=r"$\tau_p = 25.0 \mu s$")
plt.semilogy(bin_edge_log_u1_single[:-1], hist_int_log_u1_single/(u1_single.size)/(bin_edge_log_u1_single[2]-bin_edge_log_u1_single[1]),'g-', label=r"$\tau_p = 50.5 \mu s$")
plt.semilogy(bin_edge_log_u2_single[:-1], hist_int_log_u2_single/(u2_single.size)/(bin_edge_log_u2_single[2]-bin_edge_log_u2_single[1]),'r-', label=r"$\tau_p = 101.0 \mu s$")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)


############### Diff pulse width dual focused ###################



path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_Set_2_normal/Up_down_Chirp_no_RCM_CR_5050_n_az_2_13_20200830_180657_200_lines_blanked"
u1=np.load(path_temp+'/u2_single_focused.npy')
u3=np.load(path_temp+'/u2_dual_focused.npy')

path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_f_s_110MHz_tau_101us_n_az_8192_n_rg_32768/Up_down_Chirp_no_RCM_CR_10100_n_az_2_13_20200926_152610_300_lines_blanked"
u2=np.load(path_temp+'/u2_single_focused.npy')
u4=np.load(path_temp+'/u2_dual_focused.npy')


u1_single=u1[6144:10240,2048:6144]
u3_single=u3[6144:10240,2048:6144]
u2_single=u2[6144+8192:10240+8192,2048:6144]
u4_single=u4[6144+8192:10240+8192,2048:6144]


# region Probability density function of the log-intensity


# Compute the histograms
bin_num=200

hist_int_log_u1_single, bin_edge_log_u1_single = np.histogram(10*np.log10((np.abs(u1_single) ** 2)), bins=bin_num)
hist_int_log_u2_single, bin_edge_log_u2_single = np.histogram(10*np.log10((np.abs(u2_single) ** 2)), bins=bin_num)
hist_int_log_u3_single, bin_edge_log_u3_single = np.histogram(10*np.log10((np.abs(u3_single) ** 2)), bins=bin_num)
hist_int_log_u4_single, bin_edge_log_u4_single = np.histogram(10*np.log10((np.abs(u4_single) ** 2)), bins=bin_num)


# Display the probability density function

# City
plt.figure()
plt.semilogy(bin_edge_log_u1_single[:-1], hist_int_log_u1_single/(u1_single.size)/(bin_edge_log_u1_single[2]-bin_edge_log_u1_single[1]),'g-', label=r"$\tau_p = 50.5 \mu s$")
plt.semilogy(bin_edge_log_u3_single[:-1], hist_int_log_u3_single/(u3_single.size)/(bin_edge_log_u3_single[2]-bin_edge_log_u3_single[1]),'b-', label=r"$\tau_p = 50.5 \mu s (200 lines blanked)$")
plt.semilogy(bin_edge_log_u2_single[:-1], hist_int_log_u2_single/(u2_single.size)/(bin_edge_log_u2_single[2]-bin_edge_log_u2_single[1]),'r-', label=r"$\tau_p = 101.0 \mu s$")
plt.semilogy(bin_edge_log_u4_single[:-1], hist_int_log_u4_single/(u4_single.size)/(bin_edge_log_u4_single[2]-bin_edge_log_u4_single[1]),'c-', label=r"$\tau_p = 101.0 \mu s (300 lines blanked)$")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)


###### PDF conv we encoded #####

path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_Set_2_normal/Same_calibration/Up_down_Chirp_no_RCM_CR_5050_n_az_2_13_20201020_151336" #Put here the folder name wher the data is located and where the images are store

u=np.load(path_temp+'/u2_single_focused.npy')
#u1=np.load(path_temp+'/u2_dual_focused.npy')

path_temp="D:/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/B_100MHz_Set_2_normal/Same_calibration/Conv_no_RCM_CR_5050_n_az_2_13_20201020_145049"
u2=np.load(path_temp+'/u2_single_focused.npy')
#u3=np.load(path_temp+'/u2_dual_focused.npy')

u_single=u[6144:10240,2048:6144]
#u1_single=u1[6144:10240,2048:6144]
u2_single=u2[8195:10240,2048:6144]
#u3_single=u3[6144:10240,2048:6144]

# region Probability density function of the log-intensity


# Compute the histograms
bin_num=200

hist_int_log_u_single, bin_edge_log_u_single = np.histogram(10*np.log10((np.abs(u_single) ** 2)), bins=bin_num)
#hist_int_log_u1_single, bin_edge_log_u1_single = np.histogram(10*np.log10((np.abs(u1_single) ** 2)), bins=bin_num)
hist_int_log_u2_single, bin_edge_log_u2_single = np.histogram(10*np.log10((np.abs(u2_single) ** 2)), bins=bin_num)
#hist_int_log_u3_single, bin_edge_log_u3_single = np.histogram(10*np.log10((np.abs(u3_single) ** 2)), bins=bin_num)

# Display the probability density function

# City
plt.figure()
plt.semilogy(bin_edge_log_u2_single[:-1], hist_int_log_u2_single/(u2_single.size)/(bin_edge_log_u2_single[2]-bin_edge_log_u2_single[1]),'g-', label=r"Conventional SAR")
plt.semilogy(bin_edge_log_u_single[:-1], hist_int_log_u_single/(u_single.size)/(bin_edge_log_u_single[2]-bin_edge_log_u_single[1]),'r-', label=r"Waveform-encoded SAR")
#plt.semilogy(bin_edge_log_u1_single[:-1], hist_int_log_u1_single/(u1_single.size)/(bin_edge_log_u1_single[2]-bin_edge_log_u1_single[1]),'-',color='maroon', label=r"$\beta = 0.33 (dual focused)$")
#plt.semilogy(bin_edge_log_u3_single[:-1], hist_int_log_u3_single/(u3_single.size)/(bin_edge_log_u3_single[2]-bin_edge_log_u3_single[1]),'-',color='lime', label=r"$\beta = 0.7 (dual focused)$")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)