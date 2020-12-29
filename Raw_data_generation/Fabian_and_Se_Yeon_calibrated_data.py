import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
import os

# region Create the directory to save the figures 
# Parent Directories
parent_dir = "C:/Users/Fabian/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/Raw_data_generation"

# Leaf directory

directory = "Results for Fabian and Se-Yeon"
# Path
path = os.path.join(parent_dir, directory)
if not os.path.exists(path):
    os.makedirs(path)
path_temp = path + '/'

# endregion

# region Load the real TerraSAR-X data

data = readsav(
    'C:/Users/Fabian/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/Raw_data_generation/Scenes/compl_refl_city.sav')
complex_refl_city = data["compl_refl"]
complex_refl_city = complex_refl_city.astype(np.complex128)

data = readsav(
    'C:/Users/Fabian/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/Raw_data_generation/Scenes/compl_refl_forest.sav')
complex_refl_forest = data["compl_refl"]
complex_refl_forest = complex_refl_forest.astype(np.complex128)

data = readsav(
    'C:/Users/Fabian/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/Raw_data_generation/Scenes/compl_refl_lake.sav')
complex_refl_lake = data["compl_refl"]
complex_refl_lake = complex_refl_lake.astype(np.complex128)

data = readsav(
    'C:/Users/Fabian/Fabian/01-Studium/02-Master/04-ULM/07_MA/PulsCompression/Raw_data_generation/Scenes/compl_refl_town.sav')
complex_refl_town = data["compl_refl"]
complex_refl_town = complex_refl_town.astype(np.complex128)

# endregion

# region Calibrate TerraSAR-X databy a factor -51.83 dB

calib_corr_dB = -51.83  # dB
complex_refl_city = np.sqrt(10 ** (calib_corr_dB / 10)) * complex_refl_city
complex_refl_forest = np.sqrt(10 ** (calib_corr_dB / 10)) * complex_refl_forest
complex_refl_lake = np.sqrt(10 ** (calib_corr_dB / 10)) * complex_refl_lake
complex_refl_town = np.sqrt(10 ** (calib_corr_dB / 10)) * complex_refl_town

# endregion

# region Compute the mean intensity value and the corresponding standard deviation of the calibrated data

# Mean intensity value
mean_int_city = np.mean(np.abs(complex_refl_city) ** 2)
mean_int_forest = np.mean(np.abs(complex_refl_forest) ** 2)
mean_int_lake = np.mean(np.abs(complex_refl_lake) ** 2)
mean_int_town = np.mean(np.abs(complex_refl_town) ** 2)

# Standard deviation of the intensity 
std_int_city = np.std(np.abs(complex_refl_city) ** 2)
std_int_forest = np.std(np.abs(complex_refl_forest) ** 2)
std_int_lake = np.std(np.abs(complex_refl_lake) ** 2)
std_int_town = np.std(np.abs(complex_refl_town) ** 2)

# Print the results
print("Mean intensity value of the City is %s", mean_int_city, "with a standard deviation of %s", std_int_city)
print("Mean intensity value of the Forest is %s", mean_int_forest, "with a standard deviation of %s", std_int_forest)
print("Mean intensity value of the Lake is %s", mean_int_lake, "with a standard deviation of %s", std_int_lake)
print("Mean intensity value of the town is %s", mean_int_town, "with a standard deviation of %s", std_int_town)

# endregion

# region Probability density function of the intensity of the calibrated data

# Compute the histograms
bin_num=200

hist_int_city, bin_edge_city = np.histogram((np.abs(complex_refl_city) ** 2), bins=bin_num)
hist_int_forest, bin_edge_forest = np.histogram((np.abs(complex_refl_forest) ** 2), bins=bin_num)
hist_int_lake, bin_edge_lake = np.histogram((np.abs(complex_refl_lake) ** 2), bins=bin_num)
hist_int_town, bin_edge_town = np.histogram((np.abs(complex_refl_town) ** 2), bins=bin_num)


# Display the probability density function

# City
plt.figure()
plt.semilogy(bin_edge_city[:-1], hist_int_city/(complex_refl_city.size)/(bin_edge_city[2]-bin_edge_city[1]),'b-', label="City")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [linear]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)
plt.savefig(path_temp+"City_pdf_intensity")

# Forest
plt.figure()
plt.semilogy(bin_edge_forest[:-1], hist_int_forest/(complex_refl_forest.size)/(bin_edge_forest[2]-bin_edge_forest[1]),'r-', label="Forest")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [linear]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)
plt.savefig(path_temp+"Forest_pdf_intensity")

# Lake
plt.figure()
plt.semilogy(bin_edge_lake[:-1], hist_int_lake/(complex_refl_lake.size)/(bin_edge_lake[2]-bin_edge_lake[1]),'g-', label="Lake")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [linear]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)
plt.savefig(path_temp+"Lake_pdf_intensity")

# Town
plt.figure()
plt.semilogy(bin_edge_town[:-1], hist_int_town/(complex_refl_town.size)/(bin_edge_town[2]-bin_edge_town[1]),'m-', label="Town")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [linear]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)
plt.savefig(path_temp+"Town_pdf_intensity")

# endregion

# region Probability density function of the log-intensity of the calibrated data


# Compute the histograms
bin_num=200

hist_int_log_city, bin_edge_log_city = np.histogram(10*np.log10((np.abs(complex_refl_city) ** 2)), bins=bin_num)
hist_int_log_forest, bin_edge_log_forest = np.histogram(10*np.log10((np.abs(complex_refl_forest)) ** 2), bins=bin_num)
hist_int_log_lake, bin_edge_log_lake = np.histogram(10*np.log10((np.abs(complex_refl_lake) ** 2)), bins=bin_num)
hist_int_log_town, bin_edge_log_town = np.histogram(10*np.log10((np.abs(complex_refl_town) ** 2)), bins=bin_num)


# Display the probability density function

# City
plt.figure()
plt.semilogy(bin_edge_log_city[:-1], hist_int_log_city/(complex_refl_city.size)/(bin_edge_log_city[2]-bin_edge_log_city[1]),'b-', label="City")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)
plt.savefig(path_temp+"City_pdf_log_intensity")

# Forest
plt.figure()
plt.semilogy(bin_edge_log_forest[:-1], hist_int_log_forest/(complex_refl_forest.size)/(bin_edge_log_forest[2]-bin_edge_log_forest[1]),'r-', label="Fores")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)
plt.savefig(path_temp+"Forest_pdf_log_intensity")

# Lake
plt.figure()
plt.semilogy(bin_edge_log_lake[:-1], hist_int_log_lake/(complex_refl_lake.size)/(bin_edge_log_lake[2]-bin_edge_log_lake[1]),'g-', label="Lake")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)
plt.savefig(path_temp+"Lake_pdf_log_intensity")

# Town
plt.figure()
plt.semilogy(bin_edge_log_town[:-1], hist_int_log_town/(complex_refl_town.size)/(bin_edge_log_town[2]-bin_edge_log_town[1]),'m-', label="Town")
plt.grid(True,which='both',axis='both')
plt.legend(fontsize=22,loc="lower left")
plt.xlabel("Intensity [dB]",fontsize=24)
plt.ylabel("Probability density function",fontsize=24)
plt.tick_params(axis='both',labelsize=24,grid_alpha=1,grid_color='k')
figure = plt.gcf() # get current figure
figure.set_size_inches(12.74, 10.45)
plt.savefig(path_temp+"Town_pdf_log_intensity")



# endregion