'''Basic GRAPPA example using Shepp-Logan phantom.'''

import numpy as np
import matplotlib.pyplot as plt
from phantominator import shepp_logan

from grappa import grappa

if __name__ == '__main__':

    # Generate fake sensitivity maps: mps
    N = 128
    ncoils = 4
    xx = np.linspace(0, 1, N)
    x, y = np.meshgrid(xx, xx)
    mps = np.zeros((N, N, ncoils))
    mps[..., 0] = x**2
    mps[..., 1] = 1 - x**2
    mps[..., 2] = y**2
    mps[..., 3] = 1 - y**2

    # generate 4 coil phantom
    ph = shepp_logan(N)
    imspace = ph[..., None]*mps
    imspace = imspace.astype('complex')
    ax = (0, 1)
    kspace = 1/np.sqrt(N**2)*np.fft.fftshift(np.fft.fft2(
        np.fft.ifftshift(imspace, axes=ax), axes=ax), axes=ax)

    # crop 20x20 window from the center of k-space for calibration
    pd = 10
    ctr = int(N/2)
    calib = kspace[ctr-pd:ctr+pd, ctr-pd:ctr+pd, :].copy()

    # calibrate a kernel
    kernel_size = (5, 5)

    # undersample by a factor of 2 in both kx and ky
#    kspace[::2, 1::2, :] = 0  # [Nx, Ny, coils, ... ]
    kspace[1::2, ::2, :] = 0

    # reconstruct:
    # ---------------------------
    #     Inputs
    #     -------
    # K-space         : 2D multi-coil k-space data to reconstruct from. Make sure that the missing entries have exact zeros in them.
    # calib           : Calibration data (fully sampled k-space) - (kx, ky, coil).
    # kernel_size     : Size of the 2D GRAPPA kernel (kx, ky).
    # coil_axis       : Dimension holding coil data.  The other two dimensions should be image size: (sx, sy).
    # lamda           : Tikhonov regularization for the kernel calibration.
    # memmap          : Store data in Numpy memmaps.  Use when datasets are too large to store in memory.
    # memmap_filename : Name of memmap to store results in.  File is only saved if memmap=True.
    # silent          : Suppress messages to user.
    #
    #     Returns
    #     -------
    # res             : k-space data where missing entries have been filled in.
    # ---------------------------

    res = grappa(kspace, calib, kernel_size, coil_axis=-1, lamda=0.01, memmap=False)

    # Take a look
    res = np.abs(np.sqrt(N**2)*np.fft.fftshift(np.fft.ifft2(
        np.fft.ifftshift(res, axes=ax), axes=ax), axes=ax))
    res0 = np.zeros((2*N, 2*N))
    kk = 0
    for idx in np.ndindex((2, 2)):
        ii, jj = idx[:]
        res0[ii*N:(ii+1)*N, jj*N:(jj+1)*N] = res[..., kk]
        kk += 1
    plt.imshow(res0, cmap='gray')
    plt.show()
