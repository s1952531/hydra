import scipy.fft
import numpy as np

class FFT:

    def __init__(self, ncells, extent):
        self.nx = ncells[0]
        self.ny = ncells[1]
        self.nz = ncells[2]

        self.Lx = extent[0]
        self.Ly = extent[1]
        self.Lz = extent[2]

        self.kx = scipy.fft.fftfreq(self.nx , 1.0 / float(self.nx))
        self.ky = scipy.fft.fftfreq(self.ny , 1.0 / float(self.ny))
        self.kz = scipy.fft.fftfreq(self.nz+1 , 1.0 / float(self.nz+1))

        # wave numbers
        rkx_ = (2.0 * np.pi * self.kx) / self.Lx
        rky_ = (2.0 * np.pi * self.ky) / self.Ly
        rkz_ = (2.0 * np.pi * self.kz) / self.Lz

        self.rkx , self.rky , self.rkz = np.meshgrid(rkx_, rky_, rkz_, indexing='ij')


    def fftxyp2s(self, fp):
        return scipy.fft.fft2(x=fp, axes=(0, 1), overwrite_x=False, norm='ortho')

    def fftxys2p(self, fs):
        return scipy.fft.ifft2(x=fs, axes=(0, 1), overwrite_x=False, norm='ortho').real

    def diffx(self, fs):
        return 1j * self.rkx * fs

    def diffy(self, fs):
        return 1j * self.rky * fs
