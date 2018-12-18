# river_profile_spectra
Wavelet power spectra of one-dimensional time series. 

The software contained in this repository is designed to calculate power spectra of one-dimensional time series using wavelets. It calculates power (i.e. amplitude-squared) as a function of distance and wavelength (or time and frequency), and also distance-averaged power spectra, which is somewhat analogous to a Fourier transformation of the time series. A description of the methods used and examples are given in Roberts et al. (sub judice), Journal of Geophysical Research - Earth Surface. 

The code was designed to convert longitudinal river profiles (i.e. elevation as a function of distance, z(x)) into the distance-wavelength domain. These transformed time series are monotonic functions with evenly sampled elevations (i.e. $\delta_x = C$). Prior to transformation the time series are mirrored about the x and z axes in an attempt to minimise edge effects. Spectral bias was rectified using the approach described by Liu et al. (2007), in essence wavelet power spectra are normalised by scale. The code uses subroutines from the machine learning algorithms from Albanese et al. (2012), which can be found at http://mlpy.sourceforge.net. It also uses numpy and scipy routines to perform a basic Fourier Transform and for other basic tasks. The code has also been used to transform other (e.g. non-monotonic) time series successfully.

This repository includes [1] the source code to perform the wavelet transformation, it is written in python; [2] an example data file, which is the elevation of the Niger river extracted from the CGAIR SRTM digital elevation model down-sampled to 2 km (see Roberts et al., sub judice; http://srtm.csi.cgiar.org/srtmdata); [3] a plotting script to show results, this bash shell script uses routines from the Generic Mapping Tools (gmt; http://gmt.soest.hawaii.edu/projects/gmt) toolkit to do the plotting. 

A useful starting point for understanding the use of wavelets to calculate power spectra of one-dimensional time series is given by Torrence and Compo (1998). 

Albanese, D., Visintainer, R., Merler, S., Riccadonna, S., Jurman, G., Furlanello, C., 2012. mlpy: Machine Learning Python,  arXiv:1202.6548.

Liu, Y., Liang, X. S., Weinberg, R. H., 2007. Rectification of the Bias in the Wavelet Power Spectrum, Am. Met. Soc., doi: 10.1175/2007JTECHO511.1.

Torrence, C., Compo, G. P., 1998. A Practical Guide to Wavelet Analysis, Bull. Am. Met. Soc., 79(1), 61-78.

The code has been tested and benchmarked, but, as ever, use with caution and I suggest that you run your own tests to confirm veracity!  
If you have any issues with the code and/or comments contact: gareth.roberts@imperial.ac.uk.

Technical details
Python: 2.7.13
GMT: 5.1.1


