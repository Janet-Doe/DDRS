import numpy as np
import warnings
from typing import Optional, Sequence, Union

# Assume these functions are implemented elsewhere in your Python package
from lotka_voltera import lotka_volterra_complet, trace_multi_cycle

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
    """
    Compute and plot Lotka-Volterra dynamics for complete and multi-cycle cases.

    Parameters
    ----------
    a, b, c, d : float or array-like
        Positive model parameters; can be arrays of same length.
    x0, y0 : float or array-like
        Initial conditions; can be arrays of same length.
    tdeb, tfin : float, optional
        Start and end times for simulation. If None, determined automatically.
    seuil : float, optional
        Threshold value (default: -inf).
    K : float, optional
        Carrying capacity (default: inf).
    te, xe, ye : array-like, optional
        Experimental data (times, x-values, y-values).
    eqcomp : str, optional
        Constraint equation to apply, e.g. 'd = b'.
    chmax : int
        Maximum number of curves for legend display.
    """
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
        Tl, Xcl, Ycl, er, Tp, erTp = lotka_volterra_complet.lotka_volterra_complet(
            at[0], bt[0], ct[0], dt[0], x0t[0], y0t[0],
            tdeb, tfin, seuilt[0], Kt[0], te, xe, ye,
            return_cycles=False, return_equ=False
        )
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

        Tl, Xcl, Ycl, er, Tp, erTp, equl = lotka_volterra_complet.lotka_volterra_complet(
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
    trace_multi_cycle.trace_multi_cycle(Xc, Yc, equ, title=None, labels=ch)
    return Xc, Yc, equ
