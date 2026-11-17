"""Coordinate projection and catalog utilities.

Projected coordinates and depth are in km; stime is in seconds and migration time is in days.
"""
import numpy as np
from datetime import datetime
from pyproj import CRS, Transformer
from .McB import get_maxc, get_gft, b_est
from scipy.signal import savgol_filter

def apply_sg_filter( x, default_window=11, polyorder=1):
    if len(x) < 5:
        return x
    w = min(len(x)//2*2-1, default_window)
    if w < 3:
        return x
    return savgol_filter(x, w, polyorder)

def calculate_surface_dist(lat1, lon1, lat2, lon2):
    """
    Calculate the surface distance (in km) between points using the Haversine-like formula.
    
    Args:
        lat1 (float): Latitude of the first point.
        lon1 (float): Longitude of the first point.
        lat2 (np.array or float): Latitude(s) of the target point(s).
        lon2 (np.array or float): Longitude(s) of the target point(s).
        
    Returns:
        np.array or float: Distance in kilometers.
    """
    # Earth radius in meters
    R0 = 6367300.0
    
    phi1 = np.radians(lat1)
    phi2 = np.radians(lat2)
    delta_phi = np.radians(lat1 - lat2)
    delta_lambda = np.radians(lon1 - lon2)

    # Spherical law of cosines or similar projection used in original script
    # Note: Using the specific formula from the provided eqclustering logic
    dist = R0 * np.arccos(
        np.sin(phi1) * np.sin(phi2) +
        np.cos(phi1) * np.cos(phi2) * np.cos(delta_lambda)
    )
    
    return dist / 1000.0


def normalmix_1D(x, mu1i, mu2i, sigma1i, sigma2i):
    """
    Estimates parameters of a 2-component 1-D Gaussian mixture model 
    using the Expectation-Maximization (EM) algorithm.
    
    Args:
        x (np.array): The data sample to be analyzed (e.g., log10 of eta).
        mu1i, mu2i (float): Initial mean parameters for the two components.
        sigma1i, sigma2i (float): Initial standard deviation parameters.
        
    Returns:
        tuple: (cutoff_c, mu1, mu2, sigma1, sigma2, A, A2, p1, p2)
               where cutoff_c is the MLE boundary between the two distributions.
    """
    Ai = 0.5
    A2i = 0.5
    epsilon = 1.0
    epsilon0 = 0.001
    max_iter = 500
    it = 0

    while epsilon > epsilon0 and it < max_iter:
        # Expectation Step
        p1a = Ai * (1 / (np.sqrt(2 * np.pi * sigma1i**2)) * np.exp(-(((x - mu1i) / sigma1i)**2) / 2))
        p2a = A2i * (1 / (np.sqrt(2 * np.pi * sigma2i**2)) * np.exp(-(((x - mu2i) / sigma2i)**2) / 2))

        p1 = p1a / (p1a + p2a)
        p2 = p2a / (p1a + p2a)

        # Maximization Step
        A = np.sum(p1) / len(x)
        A2 = np.sum(p2) / len(x)

        mu1 = np.sum(p1 * x) / np.sum(p1)
        mu2 = np.sum(p2 * x) / np.sum(p2)

        sigma1 = np.sqrt(np.sum(p1 * (x - mu1)**2) / np.sum(p1))
        sigma2 = np.sqrt(np.sum(p2 * (x - mu2)**2) / np.sum(p2))

        # Check convergence
        diffs = [
            np.abs(Ai - A),
            np.abs(mu1i - mu1),
            np.abs(mu2i - mu2),
            np.abs(sigma1i - sigma1),
            np.abs(sigma2i - sigma2)
        ]
        epsilon = np.amax(diffs)

        # Update initial values for next iteration
        Ai, A2i, mu1i, mu2i, sigma1i, sigma2i = A, A2, mu1, mu2, sigma1, sigma2
        it += 1

    # Solve for the boundary c where the two Gaussians intersect
    # Formula: aa*x^2 + bb*x + cc = 0
    aa = -2 * sigma1**2 + 2 * sigma2**2
    bb = (2 * mu2 * 2 * sigma1**2) + (-2 * mu1 * 2 * sigma2**2)
    cc = (-mu2**2 * 2 * sigma1**2) + (mu1**2 * 2 * sigma2**2) - (np.log(sigma2 / sigma1) * 2 * sigma1**2 * 2 * sigma2**2)
    
    # Quadratic formula roots
    root_part = np.sqrt(bb**2 - 4 * aa * cc)
    c1 = (-bb + root_part) / (2 * aa)
    c2 = (-bb - root_part) / (2 * aa)

    # Select the root between the two means
    c = c1 if mu1 < c1 < mu2 else c2
    
    return c, mu1, mu2, sigma1, sigma2, A, A2, p1, p2


def build_children_map(parent_array):
    """Build a parent -> children index map once to avoid repeated scans."""
    parent_array = np.asarray(parent_array)
    finite_mask = np.isfinite(parent_array)
    parent_ids = parent_array[finite_mask].astype(int)
    child_ids = np.nonzero(finite_mask)[0]

    children_map = {}
    for parent, child in zip(parent_ids, child_ids):
        children_map.setdefault(parent, []).append(child)
    return children_map


def descend_all(parent_array, start_index, children_map=None, memo=None):
    """
    Recursively find all descendants of a parent node in the B-P tree.
    
    Args:
        parent_array (np.array): Array where the value at index i is the parent of event i.
        start_index (int): The index of the root event (e.g., a mainshock).
        children_map (dict, optional): Precomputed parent -> children lookup.
        memo (dict, optional): Memoization dictionary to speed up repeated lookups.
        
    Returns:
        list: Descendant indices including the start index.
    """
    if memo is None:
        memo = {}
    if children_map is None:
        children_map = build_children_map(parent_array)
        
    if start_index in memo:
        return memo[start_index]
        
    indices = [start_index]
    children = children_map.get(start_index, [])
    
    for child in children:
        child_indices = descend_all(parent_array, child, children_map, memo)
        indices.extend(child_indices)
            
    memo[start_index] = indices
    return indices


def convert_gps_to_utm(latitudes, longitudes, zone=None):
    """
    Convert GPS (Lat/Lon) to UTM coordinates (km).
    Fully database-independent (no EPSG).
    """
    latitudes = np.asarray(latitudes)
    longitudes = np.asarray(longitudes)

    if zone is None:
        lon_mean = np.mean(longitudes)
        zone = int(np.floor((lon_mean + 180) / 6) + 1)

    is_northern = np.mean(latitudes) >= 0
    hemi = "+north" if is_northern else "+south"

    transformer = Transformer.from_crs(
        "+proj=longlat +datum=WGS84 +no_defs",
        f"+proj=utm +zone={zone} {hemi} +datum=WGS84 +units=m +no_defs",
        always_xy=True
    )

    x, y = transformer.transform(longitudes, latitudes)
    return x / 1000.0, y / 1000.0, zone


def get_refined_b_value(magnitudes, mbin=0.1):
    """
    High-level wrapper to get a robust b-value for a given set of magnitudes.
    Combines MAXC and GFT strategies.
    """
    if len(magnitudes) < 20:
        return np.nan, np.nan
        
    # Get Mc using mixed strategy
    mc_maxc = get_maxc(magnitudes, mbin) + 0.1
    mc_gft, _, _ = get_gft(magnitudes, mbin)
    
    # Heuristic: use the more conservative Mc
    mc_final = max(mc_maxc, mc_gft) if not np.isnan(mc_gft) else mc_maxc
    
    # Estimate final b-value
    _, b, _, _ = b_est(magnitudes, mbin, mc_final)
    
    return round(b, 2), round(mc_final, 1)


def select_region(df, cluster, radius_km=100.0, period_years=5):
    """Select background seismicity based on cluster spatial/temporal extent."""
    median_x, median_y = cluster['x'].median(), cluster['y'].median()
    median_time = cluster['stime'].median()

    # Define boundaries using percentiles and fixed radius
    px = np.percentile(cluster['x'], [10, 90])
    py = np.percentile(cluster['y'], [10, 90])
    
    b_x_min = min(median_x - radius_km, px[0] - 0.5 * (px[1] - px[0]))
    b_x_max = max(median_x + radius_km, px[1] + 0.5 * (px[1] - px[0]))
    b_y_min = min(median_y - radius_km, py[0] - 0.5 * (py[1] - py[0]))
    b_y_max = max(median_y + radius_km, py[1] + 0.5 * (py[1] - py[0]))

    # Time window: period_years (seconds)
    t_win = period_years * 31536000
    b_t_min, b_t_max = median_time - t_win, median_time + t_win

    selected = df[
        (df['x'].between(b_x_min, b_x_max)) &
        (df['y'].between(b_y_min, b_y_max)) &
        (df['stime'].between(b_t_min, b_t_max))
    ]
    return selected
