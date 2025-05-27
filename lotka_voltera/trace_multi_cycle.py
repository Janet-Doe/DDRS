import numpy as np
import matplotlib.pyplot as plt
from typing import List, Optional, Sequence, Union

def trace_multi_cycle(
    Xc: Sequence[Sequence[float]],
    Yc: Sequence[Sequence[float]],
    equ: Optional[Union[Sequence[float], Sequence[Sequence[float]]]] = None,
    cylog: int = 0,
    ch: Optional[Sequence[str]] = None
) -> None:
    """
    Plot multiple Lotka-Volterra cycles in the phase plane.

    Parameters
    ----------
    Xc, Yc : sequences of sequences
        Each Xc[i], Yc[i] is an array-like trajectory for cycle i.
    equ : None, (2,), or (Nt,2)
        Equilibrium point(s). If shape is (2,), same equilibrium for all cycles;
        if (Nt,2), one per cycle; or None for no equilibrium.
    cylog : int, default 0
        If 1, use log-log axes; else linear.
    ch : sequence of str, optional
        Legend labels for each cycle. Length should match number of cycles.
    """
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
