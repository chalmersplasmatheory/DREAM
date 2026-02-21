# Various helpful routines

import numpy as np


def merge_dicts(old, new):
    """
    Update entries in the dictionary 'old' with the corresponding
    values in the dictionary 'new'. This function operates recursively.
    """
    d = old.copy()
    for k in new:
        if type(new[k]) == dict:
            if k not in d:
                d[k] = {}

            d[k] = merge_dicts(d[k], new[k])
        else:
            d[k] = new[k]

    return d


def safeTeXstring(s):
    mappings = {
        '#': r'\#',
        '_': r'\_'
    }

    for key, val in mappings.items():
        s = s.replace(key, val)
    
    return s


def scal(v):
    """
    Ensure that 'v' is a scalar.
    """
    return np.asarray(v).item()


def first_derivative(y, x, axis=0):
    """
    dy/dx along `axis` using 3-point, 2nd-order accurate stencils on a nonuniform grid.
    Uses one-sided 3-point stencils at the edges (also 2nd-order).
    
    Parameters
    ----------
    y : ndarray
        Array of values sampled at coordinates x along `axis`.
    x : (n,) array_like
        Coordinates corresponding to y along `axis` (must have n >= 3, distinct points).
        Can be increasing or decreasing (non-monotone is not supported).
    axis : int
        Axis of y that corresponds to x.
    """
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)

    n = x.size
    if n < 3:
        raise ValueError("Need at least 3 points for a 3-point derivative.")
    dx = np.diff(x)
    if np.any(dx == 0):
        raise ValueError("x contains repeated points (zero spacing).")
    # Optional: enforce monotonic (recommended for meaningful 'grid' coordinate)
    if not (np.all(dx > 0) or np.all(dx < 0)):
        raise ValueError("x must be strictly monotone (all increasing or all decreasing).")

    y_m = np.moveaxis(y, axis, -1)  # derivative axis -> last
    out = np.empty_like(y_m)

    # --- interior (i = 1..n-2): 3-point nonuniform central stencil ---
    hm = x[1:-1] - x[:-2]   # x_i - x_{i-1}
    hp = x[2:]   - x[1:-1]  # x_{i+1} - x_i

    c_im1 = -hp / (hm * (hm + hp))
    c_i   = (hp - hm) / (hm * hp)
    c_ip1 =  hm / (hp * (hm + hp))

    # broadcast coefficients across leading dimensions
    shp = (1,) * (y_m.ndim - 1) + (n - 2,)
    out[..., 1:-1] = (
        c_im1.reshape(shp) * y_m[..., :-2] +
        c_i.reshape(shp)   * y_m[..., 1:-1] +
        c_ip1.reshape(shp) * y_m[..., 2:]
    )

    # --- left edge (i = 0): use points 0,1,2 ---
    h0 = x[1] - x[0]
    h1 = x[2] - x[1]
    c0 = -(2*h0 + h1) / (h0 * (h0 + h1))
    c1 =  (h0 + h1)   / (h0 * h1)
    c2 = - h0         / (h1 * (h0 + h1))
    out[..., 0] = c0 * y_m[..., 0] + c1 * y_m[..., 1] + c2 * y_m[..., 2]

    # --- right edge (i = n-1): use points n-3,n-2,n-1 ---
    h0 = x[-1] - x[-2]
    h1 = x[-2] - x[-3]
    cN  =  (2*h0 + h1) / (h0 * (h0 + h1))
    cNm = -(h0 + h1)   / (h0 * h1)
    cNmm=  h0          / (h1 * (h0 + h1))
    out[..., -1] = cN * y_m[..., -1] + cNm * y_m[..., -2] + cNmm * y_m[..., -3]

    return np.moveaxis(out, -1, axis)

def second_derivative(y, x, axis=0):
    """
    d²y/dx² along `axis` using 3-point stencils on a nonuniform grid.
    Interior is 2nd-order accurate. Edge values use 3-point endpoint stencils
    (typically only 1st-order accurate at the boundary).
    
    Parameters
    ----------
    y : ndarray
        Array of values sampled at coordinates x along `axis`.
    x : (n,) array_like
        Coordinates corresponding to y along `axis` (must have n >= 3, distinct points).
        Can be increasing or decreasing (non-monotone is not supported).
    axis : int
        Axis of y that corresponds to x.
    """
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)

    n = x.size
    if n < 3:
        raise ValueError("Need at least 3 points for a 3-point second derivative.")
    dx = np.diff(x)
    if np.any(dx == 0):
        raise ValueError("x contains repeated points (zero spacing).")
    if not (np.all(dx > 0) or np.all(dx < 0)):
        raise ValueError("x must be strictly monotone (all increasing or all decreasing).")

    y_m = np.moveaxis(y, axis, -1)
    out = np.empty_like(y_m)

    # --- interior (i = 1..n-2): standard 3-point nonuniform second derivative ---
    hm = x[1:-1] - x[:-2]
    hp = x[2:]   - x[1:-1]

    c_im1 =  2.0 / (hm * (hm + hp))
    c_i   = -2.0 / (hm * hp)
    c_ip1 =  2.0 / (hp * (hm + hp))

    shp = (1,) * (y_m.ndim - 1) + (n - 2,)
    out[..., 1:-1] = (
        c_im1.reshape(shp) * y_m[..., :-2] +
        c_i.reshape(shp)   * y_m[..., 1:-1] +
        c_ip1.reshape(shp) * y_m[..., 2:]
    )

    # --- left edge (i = 0): use points 0,1,2 ---
    h0 = x[1] - x[0]
    h1 = x[2] - x[1]
    c0 =  2.0 / (h0 * (h0 + h1))
    c1 = -2.0 / (h0 * h1)
    c2 =  2.0 / (h1 * (h0 + h1))
    out[..., 0] = c0 * y_m[..., 0] + c1 * y_m[..., 1] + c2 * y_m[..., 2]

    # --- right edge (i = n-1): use points n-3,n-2,n-1 ---
    h0 = x[-1] - x[-2]
    h1 = x[-2] - x[-3]
    cN  =  2.0 / (h0 * (h0 + h1))
    cNm = -2.0 / (h0 * h1)
    cNmm=  2.0 / (h1 * (h0 + h1))
    out[..., -1] = cN * y_m[..., -1] + cNm * y_m[..., -2] + cNmm * y_m[..., -3]

    return np.moveaxis(out, -1, axis)



# ---------- Sensor computation ----------
def compute_S_grad_S_curv(
    y, x, dx, axis, y_floor
):
    """
    Compute S_grad and S_curv for y(x) on nonuniform center grids x.
    """
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    dx = np.asarray(dx, dtype=float)
    # reshape dx to be broadcastable to y
    shape = [1] * y.ndim
    shape[axis] = x.size
    dx = dx.reshape(shape)
    dx2 = dx**2

    # derivatives
    yp = first_derivative(y, x, axis=axis)
    ypp = second_derivative(y, x, axis=axis)

    # magnitude floor
    #y_floor = compute_floor(y, abs_floor=abs_floor, rel_floor=rel_floor, robust=robust_floor)
    denom = np.abs(y) + y_floor

    # directional scaled magnitudes
    S_grad = (dx * np.abs(yp)) / denom
    S_curv = (dx2 * np.abs(ypp)) / denom

    return S_grad, S_curv
