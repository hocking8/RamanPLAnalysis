This reads data from the Horiba XPlora and Horiba Labram.

Currently the fitting functions include: "lorentzian", "gaussian", and "voigt".

For Lorentzian, I fit A (amplitude), x0 (peak center), and w (FWHM). For Gaussian, I fit B (amplitude), mu0 (peak center), and s (variance). For Voigt, I fit A, x0, w, B, mu0, s, and C. C is the shape factor (or % that the system is Lorentzian).

When inputting guesses (list) into the fit_peaks function, you should include integer multiples of 3 or 7 depending on whether you're doing a Lorentzian/Gaussian or a Voigt. It can fit any number of peaks, but the fitting depends more strongly on your initial guesses as you increase the number of peaks. You can also input multiple sets of guesses by sending in a nested list, and the script will select the fit that results in the lowest RSS. This is helpful in cases where you fitting an entire map and there are multiple potential peak shapes.
