# river_profile_spectra
Wavelet power spectra of one-dimensional time series (e.g. longitudinal river profiles, z(x)). This repository contains code used to perform spectral analyses of longitudial river profiles presented in Roberts et al. (sub judice), Journal of Geophysical Research - Earth Surface. The transformed time series are monotonic functions with evenly sampled elevations (i.e. $\delta_x = C$). Prior to transformation the time series are mirrored about the x and z in an attempt to minimise edge effects.

This software is designed to calculate power spectra of one-dimensional time series (e.g. z(x)) using wavelets.
It uses subroutines from the machine learning algorithms published by Albanese et al. (2012), which can be found at http://mlpy.sourceforge.net. 

Albanese, D., Visintainer, R., Merler, S., Riccadonna, S., Jurman, G., Furlanello, C., 2012. mlpy: Machine Learning Python,  arXiv:1202.6548.

Its principal purpose is to convert longitudinal river profiles (i.e. elevation as a function of distance) into the distance-wavelength domain. A useful starting point for understanding the use of wavelets to calculate power spectra is given by Torrence and Compo (1998). 

Torrence, C., Compo, G. P., 1998. A Practical Guide to Wavelet Analysis, Bull. Am. Met. Soc., 79(1), 61--78.

Spectral bias has been corrected using the approach described by Liu et al. (2007). In essence the wavelet power spectrum is normalised by scale. 

Liu, Y., Liang, X. S., Weinberg, R. H., 2007. Rectification of the Bias in the Wavelet Power Spectrum, Am. Met. Soc., doi: 10.1175/2007JTECHO511.1.

This repository includes [1] the source code to perform the wavelet transform, it is written in python, [2] an example data file, which is the elevation of the Niger river extracted from the CGAIR SRTM digitial elevation model down-sampled to 2 km, see our paper for details, [3] a plotting script to show results, this bash shell script uses routines from generic mapping tools (gmt; http://gmt.soest.hawaii.edu/projects/gmt) to do the plotting. 

If you have any issues/comments contact: gareth.roberts@imperial.ac.uk.

