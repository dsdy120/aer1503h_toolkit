import sympy
import orb_mech_constants
import sys

x,y,z = sympy.symbols('x y z') # Cartesian coordinates
r,phi,theta_geo = sympy.symbols('r phi theta')

phi_as_xyz = sympy.atan2(x**2 + y**2, z)
r_as_xyz = sympy.sqrt(x**2 + y**2 + z**2)

z_as_rphi = r*sympy.cos(phi)
x_as_rphi = r*sympy.sin(phi)*sympy.cos(theta_geo)
y_as_rphi = r*sympy.sin(phi)*sympy.sin(theta_geo)

cos_phi_as_rz = z/r
sin_phi_as_rz = sympy.sqrt(1 - cos_phi_as_rz**2)

def kth_grav_zonal_harmonic(k):
    """
    Returns the expression for gravitational zonal harmonics for a given degree k.

    Parameters:
    k (int): Degree of the zonal harmonic.

    Returns:
    float: Gravitational zonal harmonic value.
    """
    if k < 0:
        raise ValueError("Degree k must be non-negative.")
    
    # Define the symbolic variables
    MU, R, JK = sympy.symbols('MU R JK') # Parameters

    # Expression for the kth zonal harmonic
    expr_k = (MU / r) * JK * (R / r)**k * sympy.legendre(k,sympy.cos(phi))

    return expr_k

def evaluate_kth_zonal_harmonic(k, mu=None, R=None, r=None, jk=None, phi=None):
    """
    Evaluate the kth zonal harmonic.

    Parameters:
    k (int): Degree of the zonal harmonic.
    mu (float): Gravitational parameter.
    R (float): Reference radius.
    r (float): Distance from the center of the body.
    jk (float): J2 coefficient.
    phi (float): Latitude.

    Returns:
    float: Evaluated value of the kth zonal harmonic.
    """
    if k < 0:
        raise ValueError("Degree k must be non-negative.")
    
    expr_k:sympy.exp = kth_grav_zonal_harmonic(k)

    for var in [mu, R, r, jk, phi]:
        if var is not None:
            expr_k = expr_k.subs(var)

    return expr_k

def pert_force_expr(zonal_expr):
    """
    Returns the expression for the perturbation force due to gravitational zonal harmonics.

    Parameters:
    zonal_expr (sympy expression): Expression for the zonal harmonic.
    mu (float): Gravitational parameter.
    R (float): Reference radius.
    r (float): Distance from the center of the body.
    jk (float): J2 coefficient.
    phi (float): Latitude.

    Returns:
    sympy expression: Expression for the perturbation force.
    """

    # Evaluate the zonal harmonic expression


    pdiff_r_x = x/r
    pdiff_r_y = y/r
    pdiff_r_z = z/r

    pdiff_phi_x = x*z/(r**3 * sympy.sin(phi))
    pdiff_phi_y = y*z/(r**3 * sympy.sin(phi))
    pdiff_phi_z = -sympy.sin(phi)/r

    pdiff_zonal_r = sympy.diff(zonal_expr, r)
    pdiff_zonal_phi = sympy.diff(zonal_expr, phi)

    
    px = -(pdiff_zonal_r * pdiff_r_x + pdiff_zonal_phi * pdiff_phi_x)
    py = -(pdiff_zonal_r * pdiff_r_y + pdiff_zonal_phi * pdiff_phi_y)
    pz = -(pdiff_zonal_r * pdiff_r_z + pdiff_zonal_phi * pdiff_phi_z)

    px = px\
        .subs(sympy.cos(phi),cos_phi_as_rz)\
        .subs(sympy.sin(phi),sin_phi_as_rz)\
        .simplify()
    
    py = py\
        .subs(sympy.cos(phi), cos_phi_as_rz)\
        .subs(sympy.sin(phi), sin_phi_as_rz)\
        .simplify()

    pz = pz\
        .subs(sympy.cos(phi), cos_phi_as_rz)\
        .subs(sympy.sin(phi), sin_phi_as_rz)\
        .simplify()
    
    # sympy.pprint(pdiff_r_x)
    # sympy.pprint(pdiff_r_y)
    # sympy.pprint(pdiff_r_z)
    # sympy.pprint(pdiff_phi_x)
    # sympy.pprint(pdiff_phi_y)
    # sympy.pprint(pdiff_phi_z)
    # sympy.pprint(pdiff_zonal_r)
    # sympy.pprint(pdiff_zonal_phi)

    # sympy.pprint(px)
    # sympy.pprint(py)
    # sympy.pprint(pz)

    return px,py,pz

if __name__ == "__main__":
    sympy.init_printing()

    for i in range(7):
        sympy.pprint(pert_force_expr(evaluate_kth_zonal_harmonic(i)))