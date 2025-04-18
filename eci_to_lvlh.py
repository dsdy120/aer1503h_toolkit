import sympy as sp

jk,mu,R,r,X,Y,Z,x,y,z,px,py,pz,i,Omega,omega,theta, u = sp.symbols(
    "jk mu R r X Y Z x y z px py pz i Omega omega theta u"
)

px_expr = (1/2)*(jk*mu*R**3/r**5)*(5*X/r*(7*Z**3/r**3-3*Z/r))
py_expr = (1/2)*(jk*mu*R**3/r**5)*(5*Y/r*(7*Z**3/r**3-3*Z/r))
pz_expr = (1/2)*(jk*mu*R**3/r**5)*(35*Z**4/r**4 - 30*Z**2/r**2+3)

p_ECI = sp.Matrix([[px_expr],[py_expr],[pz_expr]])

GRL_Z = sp.rot_axis3(Omega)
GRL_XP = sp.rot_axis1(i)
GRL_ZPP = sp.rot_axis3(u)

GRL = GRL_Z * GRL_XP * GRL_ZPP

LRG = sp.Inverse(GRL)

# sp.pprint(GRL)
# sp.pprint(LRG)

r_ECI = GRL * sp.Matrix([[r],[0],[0]])

sp.pprint(r_ECI)

X_ECI,Y_ECI,Z_ECI = r_ECI

p_ECI = sp.trigsimp(p_ECI.subs(X,X_ECI).subs(Y,Y_ECI).subs(Z,Z_ECI))

sp.pprint(p_ECI)

p_LVLH = LRG * p_ECI
p_LVLH = sp.trigsimp(sp.simplify(p_LVLH))

sp.pprint(p_LVLH)
sp.pprint(p_LVLH[0])