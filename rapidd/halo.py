from numpy import pi, sqrt, exp, piecewise, linspace
from scipy.special import erf


def eta_shm(vmin, vE, v0, vesc):

    '''analytic SHM eta arXiv:1509.01598'''

    K = v0**3 * pi * (sqrt(pi) * erf(vesc/v0) - 2 * (vesc/v0)* exp(-(vesc/v0)**2))

    def eta1(vmin): return ((v0**2 *pi )/(2*vE*K)) * ( -4 * vE * exp(-(vesc/v0)**2) + sqrt(pi) * v0 *( erf((vmin + vE)/v0) - erf((vmin-vE)/v0) ) )
    def eta2(vmin): return ((v0**2 *pi )/(2*vE*K)) * (  -2 * (vesc - vmin + vE) * exp(-(vesc/v0)**2) + sqrt(pi) * v0*( erf((vmin + vE)/v0) - erf((vmin-vE)/v0) )  )

    conds = [vmin < vesc - vE , (vesc - vE < vmin) & (vmin < vesc + vE)]

    return piecewise(vmin, conds, [eta1, eta2])


def gen_shm_table(path, vearth=250.5, v0=238, vesc=544.0):
    '''Generates the SHM eta function through analytic expression'''
    vend = vesc + vearth
    vminspace = linspace(0.0, vend, 100)
    eta_values = eta_shm(vminspace, vearth, v0, vesc)


    # Open the file for writing
    with open(path, "w") as f:
        # First line: number of entries and a "1"
        f.write(f"{len(vminspace)} 0\n")
        
        # Write vmin and eta values line by line
        for v, eta in zip(vminspace, eta_values):
            f.write(f"{v:.5E} {eta:.5E}\n")