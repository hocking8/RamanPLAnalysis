This reads data from the Horiba XPlora and Horiba Labram.

Currently the fitting functions include: "lorentzian", "gaussian", and "voigt".

For Lorentzian, I fit A (amplitude), x0 (peak center), and w (FWHM).
For Gaussian, I fit B (amplitude), mu0 (peak center), and s (variance).
For Voigt, I fit A, x0, w, B, mu0, s, and C. C is the shape factor (or % that the system is Lorentzian).

When inputting guesses into the fit_peaks function, you should include integer multiples of 3 or 7 depending
on whether you're doing a Lorentzian/Gaussian or a Voigt. It can fit any number of peaks, but the fitting
becomes more difficult as you increase the number of peaks.
