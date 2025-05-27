import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from scipy.integrate import solve_ivp
from typing import Optional, Callable, Tuple, Sequence, Union

def periode_lotka_volterra(
    a: float,
    b: float,
    c: float,
    d: float,
    x0: float,
    y0: float,
    N: int = 1000
) :
    
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

def rk4n(
    FunFcn: Callable[[float, np.ndarray], np.ndarray],
    tspan: Tuple[float, float],
    y0: Union[np.ndarray, Sequence[float]],
    N: int = 100
) -> Tuple[np.ndarray, np.ndarray]:
    t0, tfinal = tspan
    h = (tfinal - t0) / (N - 1)
    tout = np.linspace(t0, tfinal, N)

    y = np.atleast_1d(np.array(y0, dtype=float))  # Assure que y est un tableau float
    yout = np.zeros((N, y.size))
    yout[0, :] = y

    for k in range(1, N):
        t = tout[k-1]
        s1 = np.atleast_1d(FunFcn(t, y))
        s2 = np.atleast_1d(FunFcn(t + h/2, y + h*s1/2))
        s3 = np.atleast_1d(FunFcn(t + h/2, y + h*s2/2))
        s4 = np.atleast_1d(FunFcn(t + h, y + h*s3))
        y = y + (h/6)*(s1 + 2*s2 + 2*s3 + s4)
        yout[k, :] = y

    return tout, yout

def erreur(Xc, Yc, h, a, b, c, d, K, seuil, casrk4):
    if casrk4:
        Xdc = (-Xc[4:-1] + 8*Xc[3:-2] - 8*Xc[1:-4] + Xc[0:-5]) / (12*h)
        Ydc = (-Yc[4:-1] + 8*Yc[3:-2] - 8*Yc[1:-4] + Yc[0:-5]) / (12*h)
        xu, yu = Xc[2:-3], Yc[2:-3]
    else:
        Xdc = (-Xc[4:] + 8*Xc[3:-1] - 8*Xc[1:-3] + Xc[0:-4]) / (12*h)
        Ydc = (-Yc[4:] + 8*Yc[3:-1] - 8*Yc[1:-3] + Yc[0:-4]) / (12*h)
        xu, yu = Xc[2:-2], Yc[2:-2]

    if not np.isfinite(seuil):
        fx = xu * (a * (1 - xu / K) - b * yu)
        fy = yu * (-c + d * xu)
        return max(np.max(np.abs(fx - Xdc)), np.max(np.abs(fy - Ydc)))
    else:
        erx = max(np.abs(Xdc[(xu <= seuil)])) if np.any(xu <= seuil) else 0
        ind1 = (xu > seuil) & (yu <= seuil)
        if np.any(ind1):
            fx = xu[ind1] * a * (1 - xu[ind1] / K)
            erx = max(erx, np.max(np.abs(fx - Xdc[ind1])))
        ind2 = (xu > seuil) & (yu > seuil)
        if np.any(ind2):
            fx = xu[ind2] * (a * (1 - xu[ind2] / K) - b * yu[ind2])
            erx = max(erx, np.max(np.abs(fx - Xdc[ind2])))

        ery = max(np.abs(Ydc[(yu <= seuil)])) if np.any(yu <= seuil) else 0
        if np.any(ind1):
            fy = yu[ind1] * (-c)
            ery = max(ery, np.max(np.abs(fy - Ydc[ind1])))
        if np.any(ind2):
            fy = yu[ind2] * (-c + d * xu[ind2])
            ery = max(ery, np.max(np.abs(fy - Ydc[ind2])))

        return max(erx, ery)

def lotka_volterra_complet(a, b, c, d, x0, y0,
                           tdeb=None, tfin=None, seuil=-np.inf, K=np.inf,
                           te=None, xe=None, ye=None,
                           casrk4=1, h=1e-3,
                           xmin=None, xmax=None, nx=40,
                           ymin=None, ymax=None, ny=40,
                           chvit=0, cylog=0,
                           traceeq=0, tracetot=1,
                           Tp0=None):

    lvb = not np.isfinite(seuil) and not np.isfinite(K)
    testpart = (x0 == c/d and y0 == a/b)
    Tp, erTp = None, None
    if lvb and not testpart:
        if Tp0 is None:
            Tp, erTp = periode_lotka_volterra(a, b, c, d, x0, y0)
        else:
            Tp = Tp0
            erTp = 0

    if te is not None:
        tdeb, tfin = te[0], te[-1]
    else:
        print("la")
        print("tdeb tfin 1 : ", tdeb, tfin)
        if tdeb is None:
            print("tdeb is None")
            tdeb = 0
        if tfin is None:
            print("la 2")
            if lvb and not testpart:
                print("tdeb Tp : ", tdeb, Tp)
                tfin = tdeb + Tp
            else:
                raise ValueError("tfin non défini")

    print("tfin tdeb :", tfin, tdeb)
    N = int(np.ceil((tfin - tdeb) / h))
    T = np.linspace(tdeb, tfin, N)

    def F(t, Y):
        x, y = Y
        dx = 0 if x <= seuil else x * (a * (1 - x / K) - (b * y if y > seuil else 0))
        dy = 0 if y <= seuil else y * (-c + (d * x if x > seuil else 0))
        return [dx, dy]

    Y0 = [x0, y0]
    if lvb and testpart:
        T = np.array([tdeb, tfin])
        Xc = np.array([x0, x0])
        Yc = np.array([y0, y0])
        er = 0
    else:
        if casrk4:
            T, Z = rk4n(F, [T[0], T[-1]], Y0, N-1)
        else:
            sol = solve_ivp(F, [T[0], T[-1]], Y0, t_eval=T)
            T, Z = sol.t, sol.y.T
        Xc, Yc = Z[:, 0], Z[:, 1]
        er = erreur(Xc, Yc, h, a, b, c, d, K, seuil, casrk4)

    equ = None
    if traceeq or (tracetot and lvb):
        A = np.array([[a/K, b], [d, 0]])
        B = np.array([a, c])
        equ = np.linalg.solve(A, B)

    if tracetot:
        plt.figure()
        if chvit:
            if xmin is None: xmin = 0
            if xmax is None: xmax = max(Xc)
            if ymin is None: ymin = 0
            if ymax is None: ymax = max(Yc)
            x_vals = np.linspace(xmin, xmax, nx)
            y_vals = np.linspace(ymin, ymax, ny)
            X, Y = np.meshgrid(x_vals, y_vals)
            U = (X > seuil) * X * (a * (1 - X / K) - b * (Y > seuil) * Y)
            V = (Y > seuil) * Y * (-c + d * (X > seuil) * X)
            plt.quiver(X, Y, U, V, color='gray')

        if cylog:
            plt.loglog(Xc, Yc)
        else:
            plt.plot(Xc, Yc)
        if traceeq and equ is not None:
            plt.plot(*equ, 'ro')
        plt.xlabel("Proies")
        plt.ylabel("Prédateurs")
        plt.title("Trajectoire dans le plan de phase")
        plt.grid(True)

        plt.figure()
        plt.plot(T, Xc, label="Proies")
        plt.plot(T, Yc, label="Prédateurs")
        if traceeq and equ is not None:
            plt.axhline(y=equ[0], linestyle='--', color='red')
            plt.axhline(y=equ[1], linestyle='--', color='blue')
        if te is not None and xe is not None and ye is not None:
            plt.plot(te, xe, 'ro', label="Proies (mesurées)")
            plt.plot(te, ye, 'bo', label="Prédateurs (mesurés)")
        plt.xlabel("Temps")
        plt.ylabel("Effectifs")
        plt.legend()
        plt.title("Évolution temporelle")
        plt.grid(True)
        plt.show()

    return T, Xc, Yc, er, Tp, erTp, equ


def trace_multi_cycle(
    Xc: Sequence[Sequence[float]],
    Yc: Sequence[Sequence[float]],
    equ: Optional[Union[Sequence[float], Sequence[Sequence[float]]]] = None,
    cylog: int = 0,
    ch: Optional[Sequence[str]] = None
) -> None:
    
    Nt = len(Xc)
    # Convert equ to array
    equ_arr = None
    if equ is not None:
        equ_arr = np.atleast_2d(equ)
        q = equ_arr.shape[0]
    else:
        q = 0

    colors = plt.cm.hsv(np.linspace(0, 1, Nt))
    plt.figure()
    # Plot cycles
    for i in range(Nt):
        xi = np.asarray(Xc[i])
        yi = np.asarray(Yc[i])
        if cylog:
            plt.loglog(xi, yi, color=colors[i])
        else:
            plt.plot(xi, yi, color=colors[i])
    # Plot equilibria
    if q == 1:
        xeq, yeq = equ_arr[0]
        if cylog:
            plt.loglog([xeq], [yeq], 'ok')
        else:
            plt.plot([xeq], [yeq], 'ok')
    elif q > 1:
        for i in range(q):
            xeq, yeq = equ_arr[i]
            if cylog:
                plt.loglog([xeq], [yeq], 'o', color=colors[i])
            else:
                plt.plot([xeq], [yeq], 'o', color=colors[i])
    # Plot starting points
    for i in range(Nt):
        xi = np.asarray(Xc[i])
        yi = np.asarray(Yc[i])
        x0, y0 = xi[0], yi[0]
        if cylog:
            plt.loglog([x0], [y0], '*', color=colors[i])
        else:
            plt.plot([x0], [y0], '*', color=colors[i])
    # Legend if labels provided
    if ch:
        labels = []
        if q == 1:
            labels.append('Equilibre commun')
            labels.extend(ch)
        else:
            labels = ch
        plt.legend(labels, loc='best')

    plt.xlabel('effectifs des proies')
    plt.ylabel('effectifs des prédateurs')
    plt.tight_layout()
    plt.show()


def lotka_volterra_bis(
    a: Union[float, Sequence[float]],
    b: Union[float, Sequence[float]],
    c: Union[float, Sequence[float]],
    d: Union[float, Sequence[float]],
    x0: Union[float, Sequence[float]],
    y0: Union[float, Sequence[float]],
    tdeb: Optional[float] = None,
    tfin: Optional[float] = None,
    seuil: Optional[float] = None,
    K: Optional[float] = None,
    te: Optional[np.ndarray] = None,
    xe: Optional[np.ndarray] = None,
    ye: Optional[np.ndarray] = None,
    eqcomp: Optional[str] = None,
    chmax: int = 15
):
    print("tdeb : ", tdeb)
    
    # Handle defaults
    if seuil is None:
        seuil = -np.inf
    if K is None:
        K = np.inf

    # Convert all to numpy arrays
    params = [np.atleast_1d(v) for v in (a, b, c, d, x0, y0, seuil, K)]
    lengths = [v.size for v in params]
    names = ['a', 'b', 'c', 'd', 'x0', 'y0', 'seuil', 'K']

    # Determine number of parameter sets
    idx_var = [i for i, n in enumerate(lengths) if n > 1]
    if idx_var:
        p = lengths[idx_var[0]]
        # ensure all varying parameters have same length
        if any(lengths[i] != p for i in idx_var):
            raise ValueError('Multiple parameters vary with differing sizes')
    else:
        p = 1

    # Broadcast or repeat singleton parameters
    full = []
    for v in params:
        if v.size == 1:
            full.append(np.full(p, v.item()))
        else:
            full.append(v)
    at, bt, ct, dt, x0t, y0t, seuilt, Kt = full

    # Stack for legend formation
    Xt = np.vstack(full)

    # Single case
    if p == 1:
        print("LAAAAAA")
        Tl, Xcl, Ycl, er, Tp, erTp = lotka_volterra_complet(
            at[0], bt[0], ct[0], dt[0], x0t[0], y0t[0],
            tdeb, tfin, seuilt[0], Kt[0], te, xe, ye
        )
        print("LAAAAAAA 2")
        print('erreur schéma:', er)
        if not np.isnan(erTp):
            print('éventuelle erreur période:', erTp)
        return Tl, Xcl, Ycl

    # Multi-case
    Xc, Yc = [], []
    equ = np.zeros((p, 2))
    for q in range(p):
        print(f'Calcul de la courbe numéro {q+1} en cours ...')
        ai, bi, ci, di = at[q], bt[q], ct[q], dt[q]
        xi, yi = x0t[q], y0t[q]
        si, Ki = seuilt[q], Kt[q]

        # Apply constraint if any
        if eqcomp:
            # e.g. 'd = b' will assign di = bi
            local = {'a': ai, 'b': bi, 'c': ci, 'd': di,
                     'x0': xi, 'y0': yi, 'seuil': si, 'K': Ki}
            exec(eqcomp, {}, local)
            ai, bi, ci, di = local['a'], local['b'], local['c'], local['d']
            xi, yi, si, Ki = local['x0'], local['y0'], local['seuil'], local['K']

        Tl, Xcl, Ycl, er, Tp, erTp, equl = lotka_volterra_complet(
            ai, bi, ci, di, xi, yi,
            tdeb, tfin, si, Ki, te, xe, ye,
            return_cycles=True, return_equ=True
        )
        print('erreur schéma:', er)
        if not np.isnan(erTp):
            print('éventuelle erreur période:', erTp)
        Xc.append(Xcl)
        Yc.append(Ycl)
        equ[q, :] = equl

    # Consolidate equilibrium solutions
    if np.allclose(equ, equ[0, :]):
        equ = equ[0, :]

    # Legend labels
    ch = []
    if p <= chmax:
        for q in range(p):
            if len(idx_var) == 1:
                ch.append(f"{names[idx_var[0]]}={Xt[idx_var[0], q]}")
            else:
                parts = [f"{names[i]}={Xt[i, q]}" for i in idx_var]
                ch.append(', '.join(parts))

    # Plot multi-cycle
    trace_multi_cycle(Xc, Yc, equ, title=None, labels=ch)
    return Xc, Yc, equ
