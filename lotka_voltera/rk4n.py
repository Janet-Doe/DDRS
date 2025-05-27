import numpy as np
from typing import Callable, Tuple, Sequence, Union

def rk4n(
    FunFcn: Callable[[float, np.ndarray], np.ndarray],
    tspan: Tuple[float, float],
    y0: Union[np.ndarray, Sequence[float]],
    N: int = 100
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Integrate ODE system using classical 4th-order Runge-Kutta (RK4) method.

    Parameters
    ----------
    FunFcn : function
        Right-hand side function: FunFcn(t, y) returns dy/dt as array.
    tspan : tuple (t0, tfinal)
        Integration interval.
    y0 : array-like
        Initial condition vector.
    N : int, optional
        Number of time points (default: 100).

    Returns
    -------
    tout : ndarray, shape (N,)
        Time points from t0 to tfinal inclusive.
    yout : ndarray, shape (N, len(y0))
        Solution array; each row corresponds to state at tout[k].
    """
    t0, tfinal = tspan
    # Time step
    h = (tfinal - t0) / (N - 1)
    # Initialize time array
    tout = np.linspace(t0, tfinal, N)
    # Initialize solution array
    y0 = np.atleast_1d(y0)
    y = y0.astype(float)
    yout = np.zeros((N, y.size))
    yout[0, :] = y

    # RK4 loop
    for k in range(1, N):
        t = tout[k-1]
        s1 = FunFcn(t, y)
        s2 = FunFcn(t + h/2, y + h*s1/2)
        s3 = FunFcn(t + h/2, y + h*s2/2)
        s4 = FunFcn(t + h, y + h*s3)
        y = y + (h/6)*(s1 + 2*s2 + 2*s3 + s4)
        yout[k, :] = y

    return tout, yout
