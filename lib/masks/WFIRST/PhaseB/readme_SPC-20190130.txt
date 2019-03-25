----------------------------------------------------------------------
* Description of and Mask Files for WFIRST CGI Coronagraph Design "SPC-20190130" *
----------------------------------------------------------------------
Author: A.J. Riggs (Jet Propulsion Laboratory, California Institute of Technology)
Written by A.J. Riggs on January 30, 2019.

Copyright 2019 California Institute of Technology. Government sponsorship acknowledged.

The decision to implement the WFIRST mission will not be finalized until NASA’s completion of the National Environmental Policy Act (NEPA) process. This document is being made available for information purposes only.

----------------------------------------------------------------------
This file package contains :

readme_SPC-20190130.txt (this file)
pupil_SPC-20190130.fits
SPM_SPC-20190130.fits
FPM_res100_SPC-20190130.fits
LS_SPC-20190130.fits
 
----------------------------------------------------------------------
Specifications of the masks

- Design is an SPLC (shaped pupil + opaque annular focal plane mask + Lyot stop) for the WFIRST CGI.
  - Spectroscopy mode 
  - 2x65-degree field of view
  - 18% spectral bandwidth
  - 730nm center wavelength of bandpass (although the design is the same for any center wavelength since the amplitude masks are achromatic. the FPM would just have to be scaled up or down for larger or smaller wavelengths)

- TCA exit pupil (=CGI input pupil) file (pupil_SPC-20190130.fits)
  - Pixel-centered (FFT convention) in 1002x1002 array. Beam diameter is 1000 pixel widths.
  - Generated with PROPER (rectangles and circles) using a fit by A.J. Riggs to the (early) Phase B pupil "CGI 180718" from GSFC.
  - Anti-aliased (“gray”) edges.

- Shaped pupil apodizer file (SPM_SPC-20190130.fits)
  - Pixel-centered (FFT convention) in 1002x1002 array. Beam diameter is 1000 pixel widths. (Same resolution as pupil file.)
  - Binary (0 or 1) values. Non-binary amplitude values were rounded after the optimization with negligible loss in performance.
  - Starting input pupil padded (eroded) by these amounts:
    -Symmetrized obscurations about the vertical axis (as defined by stored data in the file)
    - Bulk padding:
      - +/- 0.1% of pupil diameter (D) normal to each pupil obstruction (as an example the struts are 0.4% D thicker)
      - +/- 4.0 milliradians of clocking
      - +/- 0.1% D pupil magnification
    - Additional outer diameter (OD) padding: -0.5% in pupil radius (= 1.0% D total loss in overall diameter) to block the primary mirror edge rolloff.
    - Additional central obscuration (COBS) and COBS tabs' padding: +0.254% in pupil radius (= 6.0mm/2.3631m)
    - Additional strut width padding: +/-0.123%  (= 2.9mm/2.3631m)
    - Thresholded at an amplitude of 0.99 to be binary. 

- Focal plane mask (FPM) file (FPM_res100_SPC-20190130.fits)
  - Resolution = 100 pixels per lambda_central/D. To allow user to downsample to any desired resolution.
  - Pixel-centered on 1801x1801 array. Beam diameter is 1800 pixel widths.
  - "Bowtie" occulting mask (metal on glass or through-hole Si wafer, transmissive)
  - Opening angle(s) = 2x65 degrees
  - Inner Radius = 2.60 lambda_central/D
  - Outer Radius = 9.00   lambda_central/D
  - Radius of Curvature = 0.25 lambda_central/D applied at each "corner" of the bowtie openings

- Lyot stop mask file (LS_SPC-20190130.fits)
  - Pixel-centered (FFT convention) in 1002x1002 array. Beam diameter is 1000 pixel widths. (Same resolution as pupil file.)
  - Sideways "bowtie" amplitude mask
  - Opening angle(s) = 2x90 degrees
  - Inner Radius = 38% of telescope diameter
  - Inner Radius = 92% of telescope diameter
  - No radius of curvature included in design; however, a radius of curvature of 1% diameter at "corners" does not degrade performance at all if needed for manufacturing.

- Dark Hole
  - same size of FPM opening
  - IWA = ~3.0 lambda_central/D
  - OWA = ~8.7 lambda_central/D
  - 2x65 degrees

----------------------------------------------------------------------
Ideal (No-Aberration) Performance:
  - (Main-lobe, FWHM) Core throughput = 4.05% (calculated within the half-max isophote for the main PSF lobe for a source 6 lambda/D off axis.
  - Contrast is much worse from 19-20 lambda_central/D because of the throughput drop off.
    - ~2-3e-9 raw contrast for ideal design

