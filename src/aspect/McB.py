# Legacy magnitude of completeness and b-value estimators.
"""
Author: Xin Cui
Description: Tools for calculating Magnitude of Completeness (Mc) and b-values.
Methods: MAXC (Maximum Curvature), GFT (Goodness-of-Fit), MBS (b-value Stability).
Ref: Modified from https://github.com/sachalapins/bvalues
"""

import math
import numpy as np

def calculate_fmd(mag, mbin):
    """
    Calculate the Frequency-Magnitude Distribution (FMD).
    
    Returns:
        mi (np.array): Magnitude bins.
        nbmag (np.array): Number of events in each bin.
        cumnbmag (np.array): Cumulative number of events (>= mag).
    """
    minmag = math.floor(min(mag/mbin)) * mbin
    maxmag = math.ceil(max(mag/mbin)) * mbin
    mi = np.arange(minmag, maxmag + mbin, mbin)
    nbm = len(mi)
    
    cumnbmag = np.zeros(nbm)
    for i in range(nbm):
        cumnbmag[i] = np.sum(mag > mi[i] - mbin/2)
        
    nbmag = np.abs(np.diff(np.append(cumnbmag, 0)))
    return mi, nbmag, cumnbmag

def b_est(mag, mbin, mc):
    """
    Estimate b-value using Maximum Likelihood Method (Aki, 1965).
    Includes Shi & Bolt (1982) uncertainty estimate.
    """
    mask = mag > round(mc, 1) - mbin/2
    mag_above_mc = mag[mask]
    n = len(mag_above_mc)
    
    if n < 2:
        return np.nan, np.nan, np.nan, np.nan
        
    mbar = np.mean(mag_above_mc)
    # Aki formula for b-value
    b = math.log10(math.exp(1)) / (mbar - (mc - mbin/2))
    # a-value
    a = math.log10(n) + b * mc
    # Uncertainty
    aki_unc = b / math.sqrt(n)
    shibolt_unc = 2.3 * b**2 * math.sqrt(np.sum((mag_above_mc - mbar)**2) / (n * (n-1)))
    
    return a, b, aki_unc, shibolt_unc

def get_maxc(mag, mbin):
    """Method 1: Maximum Curvature (MAXC)."""
    if len(mag) == 0: return np.nan
    mi, nbmag, _ = calculate_fmd(mag, mbin)
    maxc = mi[np.argmax(nbmag)]
    return round(maxc, 1)

def get_gft(mag, mbin, conf_levels=[90, 95]):
    """Method 2: Goodness-of-Fit Test (GFT)."""
    mi, nbmag, cumnbmag = calculate_fmd(mag, mbin)
    
    a_vals = np.zeros(len(mi))
    b_vals = np.zeros(len(mi))
    R_vals = np.zeros(len(mi))

    for i in range(len(mi)):
        m_cutoff = round(mi[i], 1)
        a, b, _, _ = b_est(mag, mbin, m_cutoff)
        if np.isnan(a):
            R_vals[i] = 100
            continue
            
        synthetic_gr = 10**(a - b * mi)
        Bi = cumnbmag[i:]
        Si = synthetic_gr[i:]
        # Eq 6: Residual calculation
        R_vals[i] = (np.sum(np.abs(Bi - Si)) / np.sum(Bi)) * 100

    mc = np.nan
    for level in conf_levels:
        fit_indices = np.where(R_vals <= (100 - level))[0]
        if len(fit_indices) > 0:
            mc = round(mi[fit_indices[0]], 1)
            break # Return the first (lowest) Mc satisfying the confidence level
            
    return mc, mi, R_vals

def get_mbs(mag, mbin, dM=0.4, min_mc=-3.0):
    """Method 3: b-value Stability (MBS)."""
    mi, nbmag, cumnbmag = calculate_fmd(mag, mbin)
    n_bins = len(mi)
    
    b_vals = np.full(n_bins, np.nan)
    shibolt_unc = np.full(n_bins, np.nan)

    for i in range(n_bins):
        if cumnbmag[i] > 1:
            _, b_vals[i], _, shibolt_unc[i] = b_est(mag, mbin, mi[i])

    # Rolling average stability check
    step_bins = round(dM / mbin)
    b_avg = np.full(n_bins, np.nan)
    is_stable = [False] * n_bins
    
    for i in range(n_bins - step_bins - 1):
        window = b_vals[i : i + step_bins + 1]
        if not np.any(np.isnan(window)):
            b_avg[i] = np.mean(window)
            is_stable[i] = abs(b_avg[i] - b_vals[i]) <= shibolt_unc[i]

    stable_points = mi[np.array(is_stable)]
    valid_stable = stable_points[stable_points > min_mc]
    
    mc = round(min(valid_stable), 1) if len(valid_stable) > 0 else np.nan
    return mc, mi, b_vals, b_avg, shibolt_unc