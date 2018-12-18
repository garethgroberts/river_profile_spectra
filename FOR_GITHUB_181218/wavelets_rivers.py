#!/usr/bin/env python

# This code is designed to calculate the power spectra of longitudinal river profiles,
# i.e. elevation as a function of distance, z(x). It does so using a continuous wavelet 
# transform based on the methods summarised by Torrence & Tompo (1998) and Liu et al. (2007).
# It makes use of the wavelet routines in the machine learning python package, mlpy. 
# The methods used and their application to geomorphic problems are discussed in our paper:
# 
# Roberts, G. G., White, N., Lodhia, B., The Generation and Scaling of Longitudinal River Profiles,
# Journal of Geophysical Research - Earth Surface, Sub Judice. 
#
# The code has been tested and benchmarked. However I suggest that you regard this software as 
# developmental and test it yourself accordingly. 
# 
# Additional useful references: 
#
# Albanese, D., Visintainer, R., Merler, S., Riccadonna, S., Jurman, G., Furlanello, C., 2012. mlpy: Machine Learning Python,  arXiv:1202.6548.
#
# Liu, Y., Liang, X. S., Weinberg, R. H., 2007. Rectification of the Bias in the Wavelet Power Spectrum, Am. Met. Soc., doi: 10.1175/2007JTECHO511.1.
#
# Torrence, C., Compo, G. P., 1998. A Practical Guide to Wavelet Analysis, Bull. Am. Met. Soc., 79(1), 61-78.
#
# Gareth Roberts, 2017-2018, UK. 
# gareth.roberts@imperial.ac.uk

from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import mlpy.wavelet as wave
from matplotlib import gridspec
from scipy.special import erf

pi=np.pi

# set time step
dt=2000.

# set wavelet function
#wf = 'dog' 
wf = 'morlet'
#wf = 'paul'

# set morlet frequency (p = omega_0, e.g. 6), paul order (p = m, e.g. 2) or dog derivative (p = m, e.g. 2)
p = 6
dj=0.1

if wf == ('morlet') :
 if p != 6:
  print 'ERROR P NOT EQUAL TO 6, ICWT WILL LOOK INCORRECT. CONTINUING ANYWAY...'
 recon_factor = (dj*(dt**0.5))/(0.776*(pi**-0.25))  # morlet: equation 11 in Torrance and Compo
 
if wf == ('dog') :
 if p == 2 :
  recon_factor = (dj*(dt**0.5))/(3.541*0.867)  # dog: equation 11 in Torrance and Compo
 if p == 6:
  recon_factor = (dj*(dt**0.5))/(1.966*0.884)
 if (p != 2) and (p !=6) :
  print 'ERROR P NOT EQUAL TO 2 or 6, ICWT WILL LOOK INCORRECT. CONTINUING ANYWAY...'
  recon_factor = (dj*(dt**0.5))/(3.541*0.867)
  print 'go'

# load data and mirror
x = np.loadtxt('one.z')
xmax = np.amax(x)
x_rev = x[::-1]
x_flip = (x * -1.) + (2. * xmax)
x_rev_flip = (x_rev * -1.) + (2. * xmax)

#x = np.concatenate((x_rev,x,x_rev_flip,x_flip,x_rev,x,x_rev_flip,x_flip),axis=0)
x = np.concatenate((x_rev,x_flip,x_rev_flip,x,x_rev,x_flip,x_rev_flip,x),axis=0)
x = x - xmax

# padding to speed up fourier transform
[x_pad, x_pad_orig] = wave.pad(x,method='zeros')
x_pad = x # no padding
xlen = len(x_pad) 
outfile = open("signal.h", 'w')
np.savetxt("signal.h", x_pad)

# fft power
fft_signal = fft.rfft(x_pad)     # Fast Fourier transform of real data
fkinv = fft.irfft(fft_signal)
fft_freq = fft.rfftfreq(xlen, d=dt)
fft_freq = fft_freq[1:(xlen/2)+1] 

fft_signal = (2. / xlen) * abs(fft_signal[1:xlen/2])
fft_power  = abs(fft_signal[1:xlen/2]) ** 2. # fft compares well to cwt
fftave = np.convolve(fft_power, np.ones(5)/5) 

# write output of fourier transform
outfile = open("fft.pow", 'w')
for i in range (0,len(fft_power)):
 print >> open('fft.pow', 'a'),  fft_freq[i], fft_power[i] , fftave[i] 
 
outfile = open("fft.inv", 'w')
for i in range (0,len(fkinv.real)):
 print >> open('fft.inv', 'a'), fkinv.real[i]
 
# calculate wavelet scales
scales = wave.autoscales(N=x_pad.shape[0], dt=dt, dj=dj, wf=wf, p=p)

# convert scales to fourier periods and output
period = wave.fourier_from_scales(scales,wf=wf,p=p)
outfile = open("scales.periods", 'w')
scale_len = len(scales)
for i in range (0, scale_len):
 print >> open("scales.periods", 'a'), i, scales[i], period[i], 1./period[i]

# perform continuous wavelet transform and output
X = wave.cwt(x=x_pad, dt=dt, scales=scales, wf=wf, p=p)
power = (np.abs(X))**2.
sumpow = power.mean(axis=1)
outfile = open("sumpower.fp", 'w')
for i in range (0, scale_len):
 print >> open('sumpower.fp', 'a'), sumpow[i] / xlen , period[i], 1./period[i], scales[i], sumpow[i] / (scales[i] * xlen )
 
# calculate cone of influence (dog wavelet)
coi = scales*(2**0.5)
outfile = open("coi.h", 'w')
np.savetxt("coi.h", coi)

# perform inverse transform to check correct scales have been used
x_icwt = wave.icwt(X=X, dt=dt, scales=scales, wf=wf,p=p)
x_icwt = x_icwt * recon_factor # equation 11 in Torrance and Compo (1998).
outfile = open("x_icwt.x", 'w')
np.savetxt("x_icwt.x", x_icwt)

# calculate mean squared error of inverse continuous wavelet transform
diff = np.sqrt((x_pad - x_icwt)**2.)
diffmean = np.mean(diff)
xmean = np.mean(abs(x_pad))
print 'Data mean = ', xmean, 'm'
print 'Transform error = ', diffmean, 'm, or', (diffmean/xmean)*100., '%'

# filtering: use only ncut largest scale (ncut longest wavelengths) coefficients for reconstruction of signal
#ncut= len(scales) / 2
ncut=87  # set, in this case, to remove all parts of signal with wavelengths < 100 km
X_numcols = len(X[0])
X_cut = X[len(scales)-ncut:len(scales),0:X_numcols]
print 'extracting coefficients at periods >', period[len(scales)-ncut], 'm (freq = ', 1./period[len(scales)-ncut], ')'
scales_cut = scales[len(scales)-ncut:len(scales)]
x_icwt_part = wave.icwt(X=X_cut, dt=dt, scales=scales_cut, wf=wf,p=p)
x_icwt_part = x_icwt_part * recon_factor # equation 11 in Torrance and Compo
outfile = open("x_icwt_part100km.x", 'w')
np.savetxt("x_icwt_part100km.x", x_icwt_part)

# calculate mean squared error of icwt_part...
diff = np.sqrt((x_pad - x_icwt_part)**2.)
diffmean = np.mean(diff)
xmean = np.mean(abs(x_pad))
print 'Data mean = ', xmean, 'm'
print 'Transform error of largest', ncut, 'scales = ', diffmean, 'm, or', (diffmean/xmean)*100., '%'

# filtering: use only ncut largest scale (ncut longest wavelengths) coefficients for # reconstruction of signal...
#ncut= len(scales) / 2
ncut2=53 # set, in this case, to remove all parts of signal with wavelengths < 1000 km
X_numcols = len(X[0])
X_cut = X[len(scales)-ncut2:len(scales),0:X_numcols]
print 'extracting coefficients at periods >', period[len(scales)-ncut2], 'm (freq = ', 1./period[len(scales)-ncut2], ')'
scales_cut = scales[len(scales)-ncut2:len(scales)]
x_icwt_part = wave.icwt(X=X_cut, dt=dt, scales=scales_cut, wf=wf,p=p)
x_icwt_part = x_icwt_part * recon_factor # equation 11 in Torrance and Compo
outfile = open("x_icwt_part1000km.x", 'w')
np.savetxt("x_icwt_part1000km.x", x_icwt_part)

# calculate mean squared error of icwt_part...
diff = np.sqrt((x_pad - x_icwt_part)**2.)
diffmean = np.mean(diff)
xmean = np.mean(abs(x_pad))
print 'Data mean = ', xmean, 'm'
print 'Transform error of largest', ncut, 'scales = ', diffmean, 'm, or', (diffmean/xmean)*100., '%'

# output power normalised by scale (e.g. Liu et al., 2007) as list of matrix elements and their values
outfile = open("file.txt", 'w')
numrows = len(power)
numcols = len(power[0])
for i in range (0,numrows):
 for j in range (0,numcols):
  print >> open('file.txt', 'a'), i, j*dt, power[i,j], period[i], 1./period[i], 1./period[len(scales)-ncut], 1./period[len(scales)-ncut2]

# remove (zero) power less than threshold (thresh)
thresh = 400.  
power[power < thresh] = 0.
sumpow = power.mean(axis=1)
outfile = open("sumpower_thresh.fp", 'w')
for i in range (0, scale_len):
 print >> open('sumpower_thresh.fp', 'a'), sumpow[i] / xlen , period[i], 1./period[i], scales[i], sumpow[i] / (scales[i] * xlen )
 
# output power normalised by scale (e.g. Liu et al., 2007) as list of matrix elements and their values
outfile = open("file_thresh.txt", 'w')
numrows = len(power)
numcols = len(power[0])
for i in range (0,numrows):
 for j in range (0,numcols):
  print >> open('file_thresh.txt', 'a'), i, j*dt, power[i,j], period[i], 1./period[i], 1./period[len(scales)-ncut]
