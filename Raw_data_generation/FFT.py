import numpy as np


def my_fft(x, savpath=None):
    # x is a vector
    x_fft = np.fft.fft(x) / len(x)

    if savpath is None:
        return x_fft
    else:
        np.save(savpath, x_fft)


def my_ifft(x, savpath=None):
    # x is a vector
    x_ifft = np.fft.ifft(x) * len(x)

    if savpath is None:
        return x_ifft
    else:
        np.save(savpath, x_ifft)


def my_fft_azimuth(x, savpath=None):
    # x is a matrix
    x_fft = np.fft.fft(x, axis=0) / x.shape[0]

    if savpath is None:
        return x_fft
    else:
        np.save(savpath, x_fft)


def my_ifft_azimuth(x, savpath=None):
    # x is a matrix
    x_ifft = np.fft.ifft(x, axis=0) * x.shape[0]

    if savpath is None:
        return x_ifft
    else:
        np.save(savpath, x_ifft)


def my_fft_range(x, savpath=None):
    # x is a matrix
    x_fft = np.fft.fft(x, axis=1) / x.shape[1]

    if savpath is None:
        return x_fft
    else:
        np.save(savpath, x_fft)


def my_ifft_range(x, savpath=None):
    # x is a matrix
    x_ifft = np.fft.ifft(x, axis=1) * x.shape[1]

    if savpath is None:
        return x_ifft
    else:
        np.save(savpath, x_ifft)
