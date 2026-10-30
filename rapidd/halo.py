from numpy import pi, sqrt, exp, piecewise, linspace
from scipy.special import erf


def eta_shm(vmin, v0, vesc, vearth):
    """
    Calculate the dark matter integrated velocity distribution eta.

    Parameters:
    vmin   : float - minimum velocity
    vesc   : float - escape velocity
    vearth : float - velocity of the Earth
    v0     : float - parameter related to the velocity distribution
    kNorm  : float - normalization constant

    Returns:
    eta    : float - integrated velocity distribution [cm^(-1) sec]
    """

    kNorm = (v0**3) * pi * (sqrt(pi) * erf(vesc / v0) - 2 * (vesc / v0) * exp(-(vesc / v0)**2))

    eta = piecewise(
        vmin,
        [
            vmin <= (vesc - vearth),
            ((vesc - vearth) < vmin) & (vmin <= (vesc + vearth)),
        ],
        [
            lambda vm: (v0**2 * pi) / (2 * vearth * kNorm)
            * ((-4) * exp(-(vesc / v0)**2) * vearth + sqrt(pi) * v0 * (erf((vm + vearth) / v0) - erf((vm - vearth) / v0))),
            lambda vm: (v0**2 * pi) / (2 * vearth * kNorm)
            * ((-2) * exp(-(vesc / v0)**2) * (vesc - vm + vearth) + sqrt(pi) * v0 * (erf(vesc / v0) - erf((vm - vearth) / v0))),
            0.0,
        ],
    )

    return eta
