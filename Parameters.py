import numpy as np
import math
from scipy.io.idl import readsav


class Parameters:

    def __init__(self):

        # %%%%%%%%%%%%%     SCENE INPUTS     %%%%%%%%%%%%%

        # Signal scene:
        # 1 -> Lake
        # 2 -> City
        # 3 -> Forest
        # 4 -> Town
        self.signal = 1

        # %%%%%%%%%%%%%     SIMULATION FLAGS     %%%%%%%%%%%%%

        # Sequence of transmitted waveforms:
        # 1 -> Cyclically-shifted chirps
        # 2 -> Up- down-chirp alternation
        self.waveform_sequence = 1

        # %%%%%%%%%%%%%     CONSTANTS     %%%%%%%%%%%%%

        # Imaginary unit
        self.j = complex(0, 1)

        # The exact value of the light speed in vacuum. [m / s]
        self.c = 2.99792458e8

        # Gravitational constant [m^3 / kg / s^2]
        self.G = 6.67384e-11

        # Earth mass [kg]
        self.M_EARTH = 5.972e24

        # Mean Earth radius [m]
        self.R_EARTH = 6371e3#6378137

        # %%%%%%%%%%%%%     PROCESSING PARAMETERS %%%%%%%%%%%%%

        self.hamming_rg = 1.  # Hamming window in range (0 -> no, 1 -> yes)
        self.hamming_rg_alpha = 0.6  # alpha coefficient for the Hamming window in range
        self.PBW = 2765  # processed Doppler bandwidth [Hz]
        self.hamming_az = 1.  # Hamming window in azimuth (0 -> no, 1 -> yes)
        self.hamming_az_alpha = 0.6  # alpha coefficient for the Hamming window in azimuth
        self.pattern_az = 1.  # azimuth antenna pattern compensation (0 -> no, 1 -> yes)
        self.plots = 1  # print the plots (0 -> no, 1 -> yes)
        self.data = 0  # print the data (0 -> no, 1 -> yes)

        # %%%%%%%%%%%%%     SYSTEM PARAMETERS     %%%%%%%%%%%%%

        # Wavelength [m]
        self.wavelength = self.c/9.65e9

        # Orbit height [m]
        self.height = 520e3

        # Chirp bandwidth [Hz]
        self.chirp_bandwidth =100e6

        # Chirp compression ratio: time - bandwidth product
        self.cr = 5050

        # Uncompressed pulse length [s]
        self.chirp_duration = self.cr / self.chirp_bandwidth

        # Range sampling frequency [Hz]
        self.rg_f_s = self.chirp_bandwidth * 1.1

        # Pulse repetition frequency [Hz]
        self.PRF = 3113

        # Antenna length [m]
        self.antenna_length = 4.8

        # Target's slant range [m]
        self.R0 = self.height+self.c/2/self.PRF#700e3

        # %%%%%%%%%%%%%     SIMULATION PARAMETERS     %%%%%%%%%%%%%

        if self.waveform_sequence == 1:
            # Sequence of short term shift-orthogonal waveforms to be used on transmit
            # Each waveform is characterized by a t_i
            self.t_i = np.zeros(2 * self.cr, dtype=np.float64)
            self.t_i = np.arange(2 * self.cr) * (np.arange(2 * self.cr) + 1) / (2 * self.chirp_bandwidth) - \
                       self.chirp_duration * (np.floor(
                (np.arange(2 * self.cr) * (np.arange(2 * self.cr) + 1) + self.cr) / (2 * self.cr)))

        elif self.waveform_sequence == 2:
            # Up- down-chirp alternation
            self.up_down_factorVec = np.array([1, -1])

        # Number of azimuth samples (typically 2^(..))
        self.n_az = 2**13#2 ** 14

        # Number of range samples --> it must be >= 2. * tau * f_s
        self.n_rg = int(2. ** np.ceil(np.log(2. * self.chirp_duration * self.rg_f_s) / np.log(2.)))#2 ** 14

        # Satellite velocity [m/s]
        self.v_s = math.sqrt(self.G * self.M_EARTH / (self.height + self.R_EARTH))

        # Ground velocity [m/s]
        self.v_g = self.v_s * self.R_EARTH / (self.height + self.R_EARTH)

        # Effective velocity [m/s]
        self.v_r = math.sqrt(self.v_s * self.v_g)

        # Azimuth time [s]
        self.t_az = (np.arange(self.n_az) - self.n_az / 2) / self.PRF

        # Range time [s] --> '2. * R0 / c0' is due to the echo delay
        self.t_rg = 2 * self.R0 / self.c + (np.arange(self.n_rg) - self.n_rg / 2) / self.rg_f_s

        # Range as a function of the azimuth time for the signal [m]
        self.R = self.R0 + self.v_r ** 2 * self.t_az ** 2 / (2 * self.R0)

        # Sinus of the azimuth angle for the signal
        self.sin_psi = (self.v_g * self.t_az / self.R0) / (1 + ((self.v_g * self.t_az / self.R0) ** 2) / 2)

        # Range as a function of the azimuth time for the nadir echo[m]
        self.h_nad = self.height + self.v_r ** 2 * self.t_az ** 2 / (2 * self.height)

        # Sinus of the azimuth angle for the signal
        self.sin_psi_nad = (self.v_g * self.t_az / self.height) / np.sqrt(1. + (self.v_g * self.t_az / self.height) ** 2.)


