"""
Name: Risa Hocking
File: RamanFitting.py
------------------
Several functions for fitting peaks in Raman. Uses the data 
structure from the RamanMap and RamanSpot classes in the
RamanReader.py file.
"""
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress
import matplotlib.pyplot as plt



def remove_linear_background(wavenumbers, intensity, ci=100):
    """
    Uses the bordering ci (int number) or so points at the end of 
    the Raman spectra to do a linear regression. Works well for scans
    where the peak takes up a lot of the spectra.

    Returns the intensity with the background substracted as well
    as the background itself.
    """
    # pull out the end points
    wn_crop = np.hstack((wavenumbers[:ci], wavenumbers[-ci:]))
    int_crop = np.hstack((intensity[:ci], intensity[-ci:]))
    
    # perform the linear regression
    m, b, r, p, se = linregress(wn_crop, int_crop)

    # compute background from this
    bg = (m * wavenumbers) + b
    return intensity - bg, bg

# this function is based on this paper: Asian Journal of Chemistry; Vol. 23, No. 12 (2011), 5229-5234
def determine_background(wavenumbers, intensity, order, eps=0.01, thresh=0.05):
    """
    Takes in the wavenumbers and intensity for a Raman peak.
    Fits the background using a polynomial of some integer order.

    wavenumbers = (N x 1) array, floats, wavenumbers from Raman
    intensity = (N x 1) array, floats, Raman intensities for one point
    order = int, polynomial order to fit
    eps = some float for fitting
    thresh = some float for fitting

    """
    # polynomical coefficient vector
    def beta(y, vdm, w):
        """
        y = intensity, (M x 1)
        vdm = Vandermonde matrix, (M x N)
        w = weight matrix, (M x M)

        Returns the polynomial coefficient vector, (N x 1). 
        """
        return np.linalg.inv(vdm.T @ w @ vdm) @ vdm.T @ w @ y
    
    # take the numerical derivative
    def numerical_deriv(wavenumbers, r):
        dwn = wavenumbers[1] - wavenumbers[0]
        dr = np.diff(r)
        d = dr/dwn
        return np.hstack((np.array([0]), d))

    # sum across some moving average window
    def sum_moving_average(i, M, d, r):
        """
        i = index of term of interest
        M = half-width of the moving window
        d = derivative
        r = residuals
        """
        # there are M points on either side of i
        if (i >= M) and ((i + M) < len(d)):
            imin = i - M
            imax = i + M
            avg = np.average(d[imin:imax])
            return np.sum((r[imin:imax] - avg)**2, axis=0)
        
        # there are only insufficient points for i + M
        elif (i >= M) and ((i + M) >= len(d)):
            imin = i - M
            avg = np.average(d[imin:])
            return np.sum((r[imin:] - avg)**2, axis=0)
        
        # there are only insufficient points for i - M
        elif (i < M) and ((i + M) < len(d)):
            imax = i + M
            avg = np.average(d[:imax])
            a = (r[:imax] - avg)**2
            return np.sum((r[:imax] - avg)**2, axis=0)
    
    # coarseness of residual
    def res_c(r, wavenumbers, M):
        # take the derivative of the residual
        d = numerical_deriv(wavenumbers, r)

        c = np.zeros(r.shape)
        for i in range(len(r)):
            sm = sum_moving_average(i, M, d, r)
            c[i] = np.sqrt(1/(2*M)) * np.sqrt(sm)

        # normalize c
        c = c / np.amax(c)
        return c
    
    # change the weight matrix
    def modify_weight(r, c, eps, thresh):
        # set the \phi(r) step function
        phir = np.where(r < 0, 1, eps).reshape((-1, 1))
        # set the \psi(c) step function
        psic = np.where(c < thresh, 1, eps).reshape(-1, 1)
        return phir @ np.reshape(psic, (1, -1))
    
    # root mean square
    def rms(x):
        return np.sqrt((1/len(x)) * np.sum(x**2))
    
    # create the vandermonde matrix for these values
    powers = np.arange(0, order+1, 1)
    pre_vdm = np.matlib.repmat(wavenumbers.reshape(len(wavenumbers), 1), 1, order+1)
    vdm = pre_vdm ** powers

    # initialize weight matrix
    w_mat = np.identity(len(wavenumbers))

    # calculate first iteration values
    b = beta(intensity, vdm, w_mat)
    res0 = intensity - (vdm @ b)
    c = res_c(res0, wavenumbers, M=10)
    w_mat = modify_weight(res0, c, eps, thresh)
    current = rms(res0)

    # set the threshold value
    good_enough = 0.1
    count = 0
    
    # iterate until rms is below the threshold
    try:
        while current > good_enough:
            count += 1
            b = beta(intensity, vdm, w_mat)
            res = intensity - (vdm @ b)
            c = res_c(res, wavenumbers, M=10)
            w_mat = modify_weight(res, c, eps, thresh)
            current = rms(res - res0)
            res0 = res
    except KeyboardInterrupt:
        print('Interrupted)')

    # obtain the estimated background
    return vdm @ b


""" Fitting profiles """
def lorentzian(x, A, x0, w):
        return A  / (((x-x0) / w)**2 + 1)
    

def gaussian(x, B, mu0, s):
    coeff = B / np.sqrt(2 * np.pi * s**2)
    exp_coeff = -(x - mu0)**2 / (2 * s**2)
    return coeff * np.exp(exp_coeff)


def voigt(x, A, x0, w, B, mu0, s, C):
    return ((1-C) * gaussian(x, B, mu0, s)) + (C * lorentzian(x, A, x0, w))


def multi_lorentzian(x, *params):
    # params should contain groups of (A, x0, w) for each Lorentzian
    n = len(params) // 3
    result = np.zeros_like(x)
    for i in range(n):
        A, x0, w = params[3*i], params[3*i+1], params[3*i+2]
        result += lorentzian(x, A, x0, w)
    return result


def multi_gaussian(x, *params):
    # params should contain groups of (B, mu0, s) for each Gaussian
    n = len(params) // 3
    result = np.zeros_like(x)
    for i in range(n):
        B, mu0, s= params[3*i], params[3*i+1], params[3*i+2]
        result += gaussian(x, B, mu0, s)
    return result


def multi_voigt(x, *params):
    # params should contain groups of (A, x0, w, B, mu0, s) for each Voigt 
    n = len(params) // 7
    result = np.zeros_like(x)
    for i in range(n):
        A, x0, w = params[7*i], params[7*i+1], params[7*i+2]
        B, mu0, s, C = params[7*i+3], params[7*i+4], params[7*i+5], params[7*i+6]
        result += voigt(x, A, x0, w, B, mu0, s, C)
    return result


def fit_peaks(wavenumbers, intensity, guesses, plot=False, method='lorentzian', fwhm_tol=0.5):
    """
    Given some data and an estimated peak position, fit
    an integer number of Voigt/Gaussian/Lorentzian curves to the data.

    Guesses is a N x M 2D nested list where M is the number of parameters for
    the fitting method you have chosen and N is the number of potential
    guesses that you want to try (if there are multiple peak shapes).

    fwhm_tol is a float that provides the tolerance value (in wavenumbers)
    for computing the full width half max of a peak.

    Can add a second guess to try in case and compare the two.
    Selects whichever guess results in a lower RSS

    Returns fit parameters, RSS, center peak, and FWHM for each peak.
    """

    # residual sum of squares to determine error
    def compute_rss(wavenumbers, fit):
        return np.sum((wavenumbers-fit)**2)
    
    # pull out the correct method for fitting
    func_dict = {'voigt': multi_voigt, 'lorentzian': multi_lorentzian, 'gaussian': multi_gaussian}
    func = func_dict[method]
    
    # initialize the fitting loop
    count = 0
    success = False # determing if this worked later
    rss = np.inf # initialize RSS as ridiculously large number
    while count < len(guesses):
        # acccount for potential that guess does not converge
        try:
            rss_t = 1E10 # some large number less than infinity
            popt_t, pcov_t = curve_fit(func, wavenumbers, intensity, p0=guesses[count], maxfev=10000)
            peak_t = func(wavenumbers, *popt_t)
            rss_t = compute_rss(intensity, peak_t)

            # preserve these fit values if they are better than previous
            if rss_t < rss:
                rss, peak, popt = rss_t, peak_t, popt_t

            success = True
        except:
            pass
    
        count += 1

    try:
        # finding peak and FWHM
        max_peak = wavenumbers[np.where(np.amax(peak) == peak)][0]
        FWHM_vec = wavenumbers[np.where(abs(np.amax(peak)/2 - peak) < fwhm_tol)]
        FWHM = abs(FWHM_vec[0] - FWHM_vec[-1])
    except:
        max_peak = np.nan
        FWHM = np.nan

    if plot and success:
        ax = plt.gca()
        # allow for the fitting to not have worked out
        # plot the fully convoluted peaks
        curve = ax.plot(wavenumbers, peak, 'k--')

        # plot the individual peaks
        if method == 'voigt':
            for i in range(len(popt) // 7):
                A, x0, w = popt[7*i], popt[7*i + 1], popt[7*i + 2]
                B, mu0, s, C = popt[7*i + 3], popt[7*i + 4], popt[7*i + 5], popt[7*i + 6]
                peak = voigt(wavenumbers, A, x0, w, B, mu0, s, C)
                curve = ax.plot(wavenumbers, peak)
                c = curve[0].get_color()
                ax.fill_between(wavenumbers, peak.min(), peak, facecolor=c, alpha=0.5)
        else:
            for i in range(len(popt) // 3):
                A, x0, w = popt[3*i], popt[3*i + 1], popt[3*i + 2]
                if method == 'lorentzian':
                    peak = lorentzian(wavenumbers, A, x0, w)
                else:
                    peak = gaussian(wavenumbers, A, x0, w)
                curve = ax.plot(wavenumbers, peak)
                c = curve[0].get_color()
                ax.fill_between(wavenumbers, peak.min(), peak, facecolor=c, alpha=0.5)

    # return the variables only if the fitting worked
    if success:
        return popt, rss, max_peak, FWHM
    else:
        return None, None, None, None


def plot_map(data, positions, cb_label=None, vmin=None, vmax=None, savefig=False, savename=None,
             allowed_range=[0, np.inf], choice = 'Raman'):
    """
    Takes in a some array of (x, y) tuples corresponding to various 
    positions, some data array with the same length as the positions 
    array, and a label for the data value in the data array.

    Allowed range is the [min, max] value for this data. Will not be
    plotted if the value is outside of this range.

    Works for all kinds of maps (rectangular, non-rectangular, and
    lines).
    """

    # collect all potential x and y values
    xs = []
    ys = []
    for x, y in positions:
        xs.append(float(x))
        ys.append(float(y))

    # remove non-unique values
    x_bare = list(set(xs))
    y_bare = list(set(ys))
    x_sort = np.asarray(sorted(x_bare))
    y_sort = np.asarray(sorted(y_bare))

    data_mat = np.zeros((len(x_sort), len(y_sort))) # initialize array to hold fit values

    # go through the data and assign everything to the correct values
    for i in range(len(positions)):
        # extract spot x and y
        x_spot = positions[i][0]
        y_spot = positions[i][1]

        # find where they are along the whole scan region
        x_ind = (np.where(x_sort == x_spot))[0][0]
        y_ind = (np.where(y_sort == y_spot))[0][0]

        data_mat[x_ind, y_ind] = data[i]

    # move the x and y min back to 0
    x_sort = x_sort - np.amin(x_sort)
    y_sort = y_sort - np.amin(y_sort)

    # remove data outside of the allowed range
    data_mat = np.where((data_mat > allowed_range[0]) & (data_mat < allowed_range[1]), data_mat, np.nan)
        
    print(f'Data ranges from {np.amin(data_mat)} to {np.amax(data_mat)}.')

    # plot a 2D grid
    fig, ax = plt.subplots(layout='constrained')
    if (len(x_sort) > 1) and (len(y_sort) > 1):
        X, Y = np.meshgrid(x_sort, y_sort)
        pc = ax.pcolormesh(y_sort, x_sort, data_mat, shading='nearest', vmin=vmin, vmax=vmax)
        plt.colorbar(pc, label=cb_label)
        ax.set_aspect('equal')
        ax.set_xlabel('X (' + r'$\mu$' + 'm)')
        ax.set_ylabel('Y (' + r'$\mu$' + 'm)')
    else:
        ax.plot(x_sort, data_mat.ravel())
        ax.set_xlabel('X (' + r'$\mu$' + 'm)')
        ax.set_ylabel(cb_label)

    if savefig:
        plt.savefig(savename, bbox_inches='tight')