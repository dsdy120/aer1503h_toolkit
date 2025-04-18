import sympy as sp

h,mu,e,h_dot,e_dot,px,py,r,theta = sp.symbols(
    "h mu e h_dot e_dot px py r theta"
)

h_dot_expr = r*py
e_dot_expr = h/mu*sp.sin(theta)*px + 1/(mu*h)*py*((h**2+mu*r)*sp.cos(theta)+mu*e*r)



rp_dot = h/(mu*(1+e))*(2/mu*h_dot - h/(1+e)*e_dot)
rp_dot = rp_dot\
    .subs(h_dot,h_dot_expr)\
    .subs(e_dot,e_dot_expr)\
    .simplify()

ra = h**2/(mu*(1+e))
ra_dot = sp.Derivative(ra,h,evaluate=True)*h_dot + sp.Derivative(ra,e,evaluate=True)*e_dot
ra_dot = ra_dot\
    .subs(h_dot,h_dot_expr)\
    .subs(e_dot,e_dot_expr)\
    .factor(py)



if __name__ == "__main__":
    sp.pprint(rp_dot)
    sp.pprint((rp_dot))

    sp.pprint(ra_dot)
    sp.pprint((ra_dot))