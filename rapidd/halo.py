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




def gen_shm_table(path, card):
    '''Generates the SHM eta function through analytic expression'''

    # --- Read this from maddm_card ---
    v0 = card['vmp'] # km/s
    vesc = card['vescape'] # km/s 

    # --- Read this from updated maddm_card --- 

    #vE_HC = [29.2, -0.1, 5.9] # Earth velocity relative to the Sun (heliocentric) at March 9 in km/s
    #vS = [11.1, 12.2, 7.3] # Solar peculiar velocity in km/s

    vE_HC = card['vEarth_mod']*[0.97986577, - 0.0033557 ,  0.19798658] # Earth velocity relative to the Sun (heliocentric) at March 9 in km/s 
    
    vS = [card['vSun_r'], card['vSun_phi'], card['cSun_theta']] # Solar peculiar velocity in km/s

    v0_vec = [0, v0, 0]

    # Calculate the module of Earth velocity relative to the galactic center
    vE_GC = sqrt(
        (vE_HC[0] + vS[0] + v0_vec[0])**2 +
        (vE_HC[1] + vS[1] + v0_vec[1])**2 +
        (vE_HC[2] + vS[2] + v0_vec[2])**2
    )

    vend = vesc + vE_GC
    vminspace = linspace(0.0, vend, 100)
    eta_values = eta_shm(vminspace, v0, vesc, vE_GC)


    # Open the file for writing
    with open(path, "w") as f:
        # First line: number of entries and a "1"
        f.write(f"{len(vminspace)} 0\n")
        
        # Write vmin and eta values line by line
        for v, eta in zip(vminspace, eta_values):
            f.write(f"{v:.5E} {eta:.5E}\n")