Code for LSD-based Analysis of Stellar Spectra 
==============================================

This package contains the code I developed for my master thesis in cooperation with dr. Andrew Tkachenko. It consists of two major parts. The first computes the least-squares deconvolution profiles of the input spectrum, using the masks given in the UserInput.dat file. The second uses the same masks and the computed LSD profiles to build a reconstruction of the observed spectrum, but with a higher SNR. Both routines are written in Fortran 90/95, although a python script is included, which provides the commands to compile & run the code, and plot the results. More information on the used algorithms can be found in:
https://ui.adsabs.harvard.edu/abs/1997MNRAS.291..658D/abstract
https://ui.adsabs.harvard.edu/abs/2010A%26A...524A...5K/abstract
https://ui.adsabs.harvard.edu/abs/2013A%26A...560A..37T/abstract

In addition, the "Data" folder contains files for a tutorial.


ACKNOWLEDGEMENTS
If the code is used for scientific publications, please cite: 
Tkachenko A., Van Reeth T., Tsymbal V., Aerts C., Kochukhov O., Debosscher J., 2013, "Denoising spectroscopic data by means of the improved least-squares deconvolution method", Astronomy & Astrophysics, Volume 560, id.A37, 14 pp.
(https://ui.adsabs.harvard.edu/abs/2013A%26A...560A..37T/abstract)


WHEN USING THE CODE FOR THE FIRST TIME:

If you want to use the overarching python script, you have to create a link to the file /Source/LSDProject.py in the main "LSDProject" folder. In GNU/Linux, this can be done by going to the LSDProject folder in the terminal and entering the command:

$ ln -s ./Source/LSDProject.py

The code can then be compiled by adding the option "-c" when running the python script:

$ python LSDProject.py -c

This has to be done twice. Error messages shown after the first compilation can be ignored. (This is due to the fact that the underlying modules have to be compiled before the executable can be compiled properly.) By default, the code uses the ifort compiler and its mkl library. If you want to use another compiler, this can be done by compiling the code "manually" or adapting the python script. The code should work for other compilers as well, though in this case, you have to (install and) use liblapack and libblas (i.e., include the options "-lblas -llapack" when compiling) instead of the mkl library.
Additional information on the use of the python script can be found by entering the command 

$ python LSDProject.py -h

in the terminal.



FOR GENERAL USE OF THE CODE:

The algorithm makes use of a large section of a normalised spectrum. In this range, there should be a large number of metal lines, but no hydrogen, helium or metal lines with damping. A good example is the section between 5000 and 5800 angstrom in the spectrum of an F-type star. The used "masks" on the other hand, are in fact line lists containing the central wavelengths, the residual fluxes at the center of the lines, and the names of the ions corresponding to the lines. (This information can for instance be obtained from the VALD website.) To allow for a more robust error estimation, there is also an option to include a file containing the variance of the observed spectrum. However, this is not necessary for the code to run. The files have to be located in the "Spectra" and "Models" subdirectories of the folder "Data". These locations then have to be filled in into the file UserInput.dat. For the tutorial this gives us:

Spectrum =      /Tutorial/synthetic_binary_070-030.dat
Variance =      /
Model =         /Tutorial/lp0000_08500_0340_0020_on_new.lin
Model =         /Tutorial/lm0003_07300_0350_0020_on_new.lin

Here we have a "/" for the variance, as we do not provide it to the code in this case. Similarly, we also have to provide the location in the "/Data/Results" folder where we want to save the results.

Result =        /Tutorial/

Other related entries include:

wavel_beg & wavel_end:   marking the edges of the studied spectral range
vel_beg & vel_end:       marking the edges of the computed profile in
                         velocity space
Regular:                 the Tikhonov regularisation parameter, which
                         allows us to "smooth out" the enhanced noise
                         resulting from the matrix inversion in the
                         algorithm. The value used in the tutorial (0.001)
                         can be considered to be "safe" in most cases. In
                         addition, this code does not smooth the LSD
                         profile itself, but the difference between the
                         LSD profile and the CCF profile of the star.
                         While this results in a slightly higher noise
                         level in the result, it does prevent the
                         smoothing of intrinsic characteristics of the
                         profile.
Nr of LSDs:              As the algorithm in this code is based on the
                         improved least-squares deconvolution (iLSD)
                         proposed by Kochukhov et al. (2010), we have the
                         possibility to compute LSD profiles for groups of
                         spectral lines with different line strengths. The
                         "separations" between the groups are provided in
                         the next line ("Line limits"). When computing a 
                         single LSD profile for our mask, the "separation
                         values" are ignored.

There is also the possibility to make use of multiple models, allowing us to compute the LSD profiles for multiple contributing stellar components, or to take into account "anomalies", such as groups of spectral lines for a pulsating star, which respond differently to the influence of pulsations.

The computation of the LSD profiles then results in the following files being created in the chosen "results" folder:

comp1.lsd:               the computed LSD profile(s). The first column
                         gives us the velocity in km/s. If N groups of
                         spectral lines were defined, based on the line
                         strength, then the following N columns give us
                         the values of the N LSD profiles in the velocity
                         domain. The final N columns give us the error
                         margins on these profiles. Naturally, if multiple
                         models were used in the computation, multiple
                         files comp<i>.lsd are created.
log:                     Contains the numerical integrals (computed with
                         the trapezium rule) of the computed LSD profiles,
                         and their sum for each individual model which was
                         used. The ratio of these integrals can be used to
                         roughly estimate the light factor of a binary, as
                         discussed in my thesis text.


***************************************************************************
***************************************************************************
***************************************************************************
    THIS PART CAN BE IGNORED FOR NOW, UNTIL I UPDATE THE CODE TO WORK WITH THE NEW GFORTRAN COMPILERS.
    
    
    The next files are also created, and form the input for the spectral reconstruction algorithm in the code.
    
    line_list.asc:           a file containing the locations of the different
                             masks which are to be used
    LSD_profile.asc:         a file containing the locations of the LSD
                             profiles which were computed
    LSC_noblocks_input.conf: file comparable to the UserInput.dat,
                             specifically for the spectral reconstruction
                             algorithm
    
    When using the python script, running the spectral reconstruction algorithm can be done by including the option "-m" or "--model":
    
    $ python LSDProject.py -m
    
    In this algorithm, the computed LSD profiles are reconvolved with the corresponding used line masks, 
    and the line strengths are optimised, so that the resulting reconstruction matches the observed spectrum 
    as well as possible. Keep in mind that, the faster the star in question is rotating, the more blending 
    there will be, and the longer this part of the code will take. (The files in the tutorial were 
    specifically chosen to illustrate the possibilities while being fairly fast. For instance HERMES spectra 
    will take a bit longer, since they have a high resolution. On the plus side, this also implies generally good results.) 
    
    If the reconstruction is succesful, the following files appear in the results folder:
    
    improved_mask.lin:       a file with the "optimised" line strengths. Keep
                             in mind that, due to a degeneracy between the
                             individual lines, these values do not have
                             physical meaning, even though the reconstructed
                             spectrum has if the algorithm was succesful!
    improved_model.dat:      the reconstruction of the observed spectrum
    O-C.dat:                 the difference between the observed and the
                             reconstructed spectrum. This file always has to
                             be checked to see if the algorithm was succesful,
                             before using the reconstructed spectrum!
    
    
    In addition, the average standard deviation for each iteration is printed to both the terminal window, and the log file. 

***************************************************************************
***************************************************************************
***************************************************************************


Finally, the computed LSD profiles and (if included) the reconstructed spectrum can be plotted by adding the option "-p" or "--plot" when running the python script.

$ python LSDProject.py -p


Questions, remarks and information on encountered bugs can always be send to timothy.vanreeth@kuleuven.be.


Good luck with your analyses!

Timothy Van Reeth
