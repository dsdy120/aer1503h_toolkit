import sympy as sp
from eci_to_lvlh import LRG,GRL,p_LVLH,px,py
from gauss_variational_eqns import e_dot_expr

sp.pprint(e_dot_expr)

sp.pprint(e_dot_expr.subs(px,p_LVLH[0]).subs(py,p_LVLH[1]).simplify())