import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Sequence, Union, Tuple

# Assume these helper functions exist in yourmodule
from lotka_voltera import rk4n, periode_lotka_volterra


def lotka_volterra_complet(
    a: float,
    b: float,
    c: float,
    d: float,
    x0: float,
    y0: float,
    tdeb: Optional[float] = None,
    tfin: Optional[float] = None,
    seuil: Optional[float] = None,
    K: Optional[float] = None,
    te: Optional[np.ndarray] = None,
    xe: Optional[np.ndarray] = None,
    ye: Optional[np.ndarray] = None,
    casrk4: int = 1,
    h: float = 1e-3,
    xmin: Optional[float] = None,
    xmax: Optional[float] = None,
    nx: int = 40,
    ymin: Optional[float] = None,
    ymax: Optional[float] = None,
    ny: int = 40,
    chvit: int = 0,
    cylog: int = 0,
    traceeq: int = 0,
    tracetot: int = 1,
    Tp0: Optional[float] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, Optional[float], Optional[float], np.ndarray]:
    """
    Compute and optionally plot trajectories of the complete Lotka-Volterra model.

    Returns:
        T: times array
        Xc, Yc: prey and predator trajectories
        er: numerical error estimate
        Tp, erTp: period and its error (for the unthresholded model)
        equ: equilibrium point [x_eq, y_eq]
    """
    # Default parameters
    if seuil is None:
        seuil = -np.inf
    if K is None:
        K = np.inf
    # Logarithmic velocity field not compatible
    if cylog == 1 and chvit == 1:
        raise ValueError("Cannot plot velocity field on log scale")

    # Base LV model flag
    lvb = not np.isfinite(seuil) and not np.isfinite(K)
    # special initial at equilibrium
    testpart = (x0 == c/d and y0 == a/b)

    # Initialize period outputs
    Tp = None
    erTp = None
    if lvb:
        if Tp0 is None and not testpart:
            Tp, erTp = periode_lotka_volterra(a, b, c, d, x0, y0)
        else:
            Tp = Tp0
            erTp = None

    # Time boundaries
    if te is not None:
        tdeb, tfin = te[0], te[-1]
    else:
        tdeb = 0.0 if tdeb is None else tdeb
        if tfin is None:
            if lvb and not testpart and Tp is not None:
                tfin = tdeb + Tp
            else:
                raise ValueError("tfin must be specified")

    # Build time array
    if not lvb:
        N = int(np.ceil((tfin - tdeb)/h))
        T = np.linspace(tdeb, tfin, N)
    else:
        if testpart:
            T = np.array([tdeb, tfin])
        else:
            N = int(np.ceil(Tp/h)) - 1
            T = np.linspace(tdeb, tdeb + Tp - h, N)

    # Velocity field for plotting
    if chvit:
        def Fx(x, y): return (x>seuil)*x*(a*(1-x/K) - b*(y>seuil)*y)
        def Fy(x, y): return (y>seuil)*y*(-c + d*(x>seuil)*x)

    # ODE right-hand side
    def F(t, Y):
        x, y = Y
        dx = (x>seuil)*x*(a*(1-x/K) - b*(y>seuil)*y)
        dy = (y>seuil)*y*(-c + d*(x>seuil)*x)
        return np.array([dx, dy])

    # Equilibrium
    equ = None
    if traceeq and tracetot:
        equ = np.linalg.solve(np.array([[a/K, b],[d, 0]]), np.array([a, c]))

    # Handle trivial equilibrium-case trajectory
    if lvb and testpart:
        Xc = np.array([x0, x0])
        Yc = np.array([y0, y0])
        er = 0.0
    else:
        # Integrate
        if casrk4:
            T, Z = rk4n(F, (T[0], T[-1]), np.array([x0, y0]), len(T))
            if te is not None:
                Zb = None  # can interpolate if needed
        else:
            from scipy.integrate import ode45  # placeholder for MATLAB ode45
            sol = ode45(F, T, np.array([x0, y0]))
            Z = sol.y.T
            if te is not None:
                # interpolate solution at te
                Zb = np.vstack([np.interp(te, sol.t, sol.y[i,:]) for i in (0,1)]).T
        Xc = Z[:,0]
        Yc = Z[:,1]

        # Error estimate
        def erreur(Xc_arr, Yc_arr):
            h0 = (T[-1] - T[0])/(len(T)-1)
            # Finite-difference derivative
            if casrk4:
                Xdc = (-Xc_arr[4:-1] + 8*Xc_arr[3:-2] - 8*Xc_arr[1:-4] + Xc_arr[0:-5])/(12*h0)
                Ydc = (-Yc_arr[4:-1] + 8*Yc_arr[3:-2] - 8*Yc_arr[1:-4] + Yc_arr[0:-5])/(12*h0)
                xu = Xc_arr[2:-3]
                yu = Yc_arr[2:-3]
            else:
                Xdc = (-Xc_arr[4:] + 8*Xc_arr[3:-1] - 8*Xc_arr[1:-3] + Xc_arr[0:-4])/(12*h0)
                Ydc = (-Yc_arr[4:] + 8*Yc_arr[3:-1] - 8*Yc_arr[1:-3] + Yc_arr[0:-4])/(12*h0)
                xu = Xc_arr[2:-2]
                yu = Yc_arr[2:-2]
            if not np.isfinite(seuil):
                errx = np.max(np.abs(xu*(a*(1-xu/K)-b*yu) - Xdc))
                erry = np.max(np.abs(yu*(-c + d*xu) - Ydc))
                return max(errx, erry)
            # piecewise
            erx = ery = 0.0
            # x errors
            inds = xu <= seuil
            if inds.any(): erx = max(erx, np.max(np.abs(Xdc[inds])))
            inds = (xu>seuil)&(yu<=seuil)
            if inds.any(): erx = max(erx, np.max(np.abs(xu[inds]*(a*(1-xu[inds]/K))-Xdc[inds])))
            inds = (xu>seuil)&(yu>seuil)
            if inds.any(): erx = max(erx, np.max(np.abs(xu[inds]*(a*(1-xu[inds]/K)-b*yu[inds]) - Xdc[inds])))
            # y errors
            inds = yu <= seuil
            if inds.any(): ery = max(ery, np.max(np.abs(Ydc[inds])))
            inds = (xu>seuil)&(yu<=seuil)
            if inds.any(): ery = max(ery, np.max(np.abs(yu[inds]*c + Ydc[inds])))
            inds = (xu>seuil)&(yu>seuil)
            if inds.any(): ery = max(ery, np.max(np.abs(yu[inds]*(-c + d*xu[inds]) - Ydc[inds])))
            return max(erx, ery)

        er = erreur(Xc, Yc)

    # Plotting
    if tracetot:
        # Phase plane
        if chvit:
            if xmin is None:
                xmin, xmax = 0, np.max(Xc)
                ymin, ymax = 0, np.max(Yc)
            xg = np.linspace(xmin, xmax, nx)
            yg = np.linspace(ymin, ymax, ny)
            Xg, Yg = np.meshgrid(xg, yg)
            U = (Xg>seuil)*Xg*(a*(1-Xg/K) - b*(Yg>seuil)*Yg)
            V = (Yg>seuil)*Yg*(-c + d*(Xg>seuil)*Xg)
            plt.figure(); plt.quiver(Xg, Yg, U, V)
        # Trajectory
        plt.figure();
        if cylog:
            plt.loglog(Xc, Yc, '-');
        else:
            plt.plot(Xc, Yc, '-');
        plt.plot(x0, y0, '*r')
        if traceeq:
            plt.plot(equ[0], equ[1], 'or')
        plt.xlabel('prey'); plt.ylabel('predator');

        # Time series
        plt.figure();
        plt.plot(T, Xc, 'r', label='prey')
        plt.plot(T, Yc, 'b', label='predator')
        if te is not None:
            plt.plot(te, xe, 'r*', te, ye, 'b*')
        if traceeq:
            plt.hlines([equ[0], equ[1]], T[0], T[-1], colors=['r','b'], linestyles='--')
        plt.xlabel('time'); plt.ylabel('population'); plt.legend()

    return T, Xc, Yc, er, Tp, erTp, equ
