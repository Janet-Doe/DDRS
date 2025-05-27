import numpy as np
from scipy.optimize import root_scalar

def periode_lotka_volterra(
    a: float,
    b: float,
    c: float,
    d: float,
    x0: float,
    y0: float,
    N: int = 1000
) :
    """
    Determine the period T of the base Lotka-Volterra model and error indicator.

    Equations (no thresholds or capacity):
        x'(t) = x * (a - b * y)
        y'(t) = y * (-c + d * x)

    Parameters
    ----------
    a, b, c, d : float
        Positive model parameters.
    x0, y0 : float
        Initial conditions (must not be the equilibrium c/d, a/b).
    N : int, optional
        Number of discretization points (default 1000).

    Returns
    -------
    T : float
        The computed period.
    er : float
        Maximum residual error of the nonlinear solves.
    """
    # Check not at equilibrium
    if x0 == c/d and y0 == a/b:
        raise ValueError("Cannot compute period starting at equilibrium")

    # Prepare initial polar coordinates
    p0 = np.log(d * x0 / c)
    q0 = np.log(b * y0 / a)
    th0, r0 = np.arctan2(p0, q0), np.hypot(p0, q0)

    # Hamiltonian-like function to satisfy
    def H(r, theta):
        return (a * (r * np.sin(theta) - np.exp(r * np.sin(theta))) +
                c * (r * np.cos(theta) - np.exp(r * np.cos(theta))) -
                a * (q0 - np.exp(q0)) -
                c * (p0 - np.exp(p0)))

    # Discretize theta
    thetas = th0 + np.linspace(0, 2 * np.pi, N)
    rs = np.empty(N)
    errs = np.empty(N)
    rs[0] = r0
    errs[0] = abs(H(r0, thetas[0]))

    # Solve for r at each theta
    for k in range(1, N):
        sol = root_scalar(lambda r: H(r, thetas[k]),
                          bracket=[0, rs[k-1] * 2],
                          x0=rs[k-1])
        rk = sol.root
        rs[k] = rk
        errs[k] = abs(H(rk, thetas[k]))

    # Compute period integral
    h = 2 * np.pi / (N - 1)
    sinth = np.sin(thetas)
    costh = np.cos(thetas)
    integrand = rs / (a * sinth * (1 - np.exp(rs * sinth)) +
                      c * costh * (1 - np.exp(rs * costh)))
    T = -h * np.trapz(integrand, dx=1)
    er = np.max(errs)
    return T, er
