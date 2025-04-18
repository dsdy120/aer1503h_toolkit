import sympy as sp
import orb_mech_constants

a,h,e,theta,pr,ps,mu,r,j2,R,x,y,z,n,omega,i = sp.symbols(
    "a h e theta pr ps mu r j2 R x y z n omega i"
)

da_dt_expr = 2*a**2/h*(e*sp.sin(theta)*pr+h**2*ps/(mu*r))
# pr_expr = (3/2)*(j2*mu*R**4/r**4)*(x/r)*(5*(z/r)**2-1)
# ps_expr = (3/2)*(j2*mu*R**4/r**4)*(y/r)*(5*(z/r)**2-1)
pr_expr = (3/2)*(j2*mu*R**2/r**4)*(1-3*(sp.sin(i))**2*(sp.sin(omega+theta))**2)
ps_expr = (3/2)*(j2*mu*R**2/r**4)*((sp.sin(i))**2*(sp.sin(2*(omega+theta))))

da_dt_expr = da_dt_expr\
    .subs(pr,pr_expr)\
    .subs(ps,ps_expr)\
    # .factor()
    # .simplify()

sp.pprint(da_dt_expr)

mean_da_dt_expr = (n/(2*sp.pi))*sp.integrate(da_dt_expr*r**2/h,(theta,0,2*sp.pi))

sp.pprint(mean_da_dt_expr)