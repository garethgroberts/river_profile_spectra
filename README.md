# river_profile_spectra
Wavelet power spectra of one-dimensional time series. 

The software contained in this repository is designed to calculate power spectra of one-dimensional time series using wavelets. It calculates power (i.e. amplitude-squared) as a function of distance and wavelength (or time and frequency), and also distance-averaged power spectra, which is somewhat analogous to a Fourier transformation of the time series. A description of the methods used and examples are given in: 

Roberts, G. G., White, N., Lodhia, B., 2019. The Generation and Scaling of Longitudinal River Profiles, Journal of Geophysical Research - Earth Surface, doi:10.1029/2018JF004796. 
https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JF004796

And relatedly, a method to calculate similarities and disparities between time series (e.g. river profiles) across the scales of interest using a cross wavelet approach is given in:

Roberts, G. G., 2019. Scales of Similarity and Disparity Between Drainage Networks, Geophysical Research Letters, doi: 10.1029/2019GL082446. 
https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019GL082446

The code was designed to convert longitudinal river profiles (i.e. elevation as a function of distance, z(x)) into the distance-wavelength domain. The time series being transformed are monotonic functions with evenly sampled elevations (i.e. $\delta_x = C$). Prior to transformation the time series are mirrored about the x and z axes in an attempt to minimise edge effects. Spectral bias was rectified using the approach described by Liu et al. (2007), in essence wavelet power spectra are normalised by scale. The code uses subroutines from the machine learning algorithms from Albanese et al. (2012), which can be found at http://mlpy.sourceforge.net, which are based on the methods described by Torrence and Compo (1998). It also uses numpy and scipy routines to perform a basic Fourier transform and for other basic computational tasks. The code has been used to transform other (e.g. non-monotonic) time series successfully (see Roberts & Mannion, 2019, Scientific Reports, https://doi.org/10.1038/s41598-019-42538-7).



This repository includes [1] the source code to perform the wavelet transformation, it is written in python; [2] an example data file, which is the elevation of the Niger river extracted from the CGAIR SRTM digital elevation model down-sampled to 2 km (one.z, see Roberts et al., sub judice; http://srtm.csi.cgiar.org/srtmdata); [3] a plotting script to show results, this bash shell script uses routines from the Generic Mapping Tools (gmt; http://gmt.soest.hawaii.edu/projects/gmt) toolkit to do the plotting. For completeness/benchmarking the plotting script is a pared down version of the one used to generate Figure 3 in Roberts et al. (2019). Note that the full resolution (~90 m) river profile is also included for completeness (obs_river), and for use in the plotting script. 

The code has been tested and benchmarked. However, I suggest that you regard it as developmental that you run your own tests to confirm veracity. If you have any issues with running the code and/or comments contact gareth.roberts@imperial.ac.uk.

Software used to generate code: Python 2.7.13, GMT 5.1.1, mlpy 3.5.0. It was developed and run on OSX 10.12.6 (macOS Sierra), with DeveloperTools installed, and should be portable to most *nix systems/python environments. 

To run the code and plot output, assuming that you have already installed the dependencies (e.g. python, numpy, scipy, mlpy, gmt) try:

> python wavelets_rivers.py

> ./plot_wave_power_rivers.gmt

If you want to run the code multiple times, say to test different mother wavelets, you can move output around using 

> ./rename_output.sh.

Examples of output for the DOG and Morlet mother wavelets are included in ./wavelet_tests. Note, however, that for compactness I've not uploaded the power spectral maps (i.e. file*.txt, or surf.grd). You will need to run the code first if you want to see that output (or email me). 

Some useful references:

Albanese, D., Visintainer, R., Merler, S., Riccadonna, S., Jurman, G., Furlanello, C., 2012. mlpy: Machine Learning Python,  arXiv:1202.6548.

Liu, Y., Liang, X. S., Weinberg, R. H., 2007. Rectification of the Bias in the Wavelet Power Spectrum, Am. Met. Soc., doi: 10.1175/2007JTECHO511.1.

Torrence, C., Compo, G. P., 1998. A Practical Guide to Wavelet Analysis, Bull. Am. Met. Soc., 79(1), 61-78.

Gareth Roberts, 2017-2019, UK. 


