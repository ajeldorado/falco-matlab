----------------------------------------------------------------------
* Description of and Mask Files for WFIRST CGI Wide-FOV Coronagraph Design "SPC-20181220" *
----------------------------------------------------------------------
Author: A.J. Riggs (Jet Propulsion Laboratory, California Institute of Technology)
Written on 2019-02-04 by A.J. Riggs.
Updated on 2019-02-05 by A.J. Riggs.

Copyright 2019 California Institute of Technology. Government sponsorship acknowledged.

The decision to implement the WFIRST mission will not be finalized until NASA’s completion of the National Environmental Policy Act (NEPA) process. This document is being made available for information purposes only.

----------------------------------------------------------------------
This file package contains :

readme_SPC-20181220.txt (this file)
pupil_SPC-20181220_1k.fits
SPM_SPC-20181220_1000_rounded9_gray.fits 
SPM_SPC-20181220_1000upsamp3x_binary.fits
FPM_res50_SPC-20181220.fits
LS_SPC-20181220_1k.fits
 
----------------------------------------------------------------------
Specifications of the masks

- Design is an SPLC (shaped pupil + opaque annular focal plane mask + Lyot stop) for the WFIRST CGI.
  - Wide field of view (FOV) imaging mode 
  - 360-degree field of view
  - 10% spectral bandwidth
  - 825nm center wavelength of bandpass (although the design is the same for any center wavelength since the amplitude masks are achromatic. the FPM would just have to be scaled up or down for larger or smaller wavelengths)

- TCA exit pupil (=CGI input pupil) file (pupil_SPC-20181220_1k.fits)
  - Pixel-centered (FFT convention) in 1002x1002 array. Beam diameter is 1000 pixel widths.
  - Generated with PROPER (rectangles and circles) using a fit by A.J. Riggs to the (early) Phase B pupil "CGI 180718" from GSFC.
  - Anti-aliased (“gray”) edges.

- Shaped pupil apodizer files 
  - File at MODELING resolution (1000x1000): SPM_SPC-20181220_1000_rounded9_gray.fits
    - Pixel-centered (FFT convention) in 1002x1002 array. Beam diameter is 1000 pixel widths. (Same resolution as pupil file.)
    - Nearly binary (0 or 1) values. Final apodizer solution was rounded to nearest 1/9. 99.46% of values are binary; the remaining 0.46% of pixels are each converted to a 3x3 sub-array with binary values.
  - File at MANUFACTURING resolution (3000x3000): SPM_SPC-20181220_1000upsamp3x_binary.fits
      - Pixel-centered (FFT convention) in 3003x3003 array. Beam diameter is 3000 pixel widths. 
      - Binary (0 or 1) values. Upsampled 3x in a smart way from the 1000x1000 SPM file to make non-binary pixels into binary 3x3 sub-arrays that avoid free-floating sub-pixels.
  - Starting input pupil padded (eroded) by these amounts:
    -Symmetrized obscurations about the vertical axis (as defined by stored data in the file)
    - Bulk padding:
      - +/- 0.2% of pupil diameter (D) normal to each pupil obstruction (as an example the struts are 0.4% D thicker)
      - +/- 4.0 milliradians of clocking
      - +/- 0.2% D pupil magnification
    - Additional outer diameter (OD) padding: -0.5% in pupil radius (= 1.0% D total loss in overall diameter) to block the primary mirror edge rolloff.
    - Additional central obscuration (COBS) and COBS tabs' padding: +0.254% in pupil radius (= 6.0mm/2.3631m)
    - Additional strut width padding: +/-0.123%  (= 2.9mm/2.3631m)
    - Thresholded at an amplitude of 0.95 to be binary. 

- Focal plane mask (FPM) file (FPM_res50_SPC-20181220.fits)
  - This is mask is generated with two concentric circles. You may want just to generate it yourself rather than downsampling the file provided here.
  - Resolution = 50 pixels per lambda_central/D. To allow user to downsample to any desired resolution.
  - Pixel-centered on 2002x2002 array. Mask opening diameter is 2000 pixel widths.
  - Annular-opening occulting mask (opaque metal on glass, transmissive)
  - Opening angle = 360 degrees
  - Inner Radius = 5.40 lambda_central/D
  - Outer Radius = 20.00   lambda_central/D

- Lyot stop mask file (LS_SPC-20181220_1k.fits)
  - Pixel-centered (FFT convention) in 1002x1002 array. Beam diameter is 1000 pixel widths. (Same resolution as pupil file.)
  - Amplitude-only mask with anti-aliased edges
  - Obscurations symmetrized about the vertical axis (as given in the file)
  - Inner Radius = 38% of telescope diameter, concentric with outer diameter
  - Inner Radius = 91% of telescope diameter
  - Strut widths = 3.2% of telescope diameter. Because the SPM's image is spatially filtered at the Lyot (pupil) plane, the struts do not have to be oversized or aligned exactly. The Lyot stop struts for SPC-20181220 are primarily there to give better structural support. 

- Dark Hole
  - Same size of FPM opening
  - IWA = ~5.7 lambda_central/D
  - OWA = ~19.7 lambda_central/D
  - 360 degrees

----------------------------------------------------------------------
Ideal (No-Aberration) Performance:
  - (FWHM) Core throughput = 4.59% (calculated within the half-max isophote for a source 13 lambda/D off axis.
  - Contrast is much worse from 19-20 lambda_central/D because of the throughput drop off.
  - ~1e-9 raw contrast from IWA to OWA for ideal design

