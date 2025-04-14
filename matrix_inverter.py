import sympy

h,e,i = sympy.symbols('h e i')

E1 = (1-e)/(1+e)
E2 = e*h/(1+e)**2
E3 = sympy.cos(i)
E4 = h*sympy.sin(i)

A = sympy.Matrix([
    [   0,  0,-E1,-E3,  0, -1],
    [   0,  0,-E2,  0,  0,  0],
    [  E1, E2,  0,  0,  0,  0],
    [  E3,  0,  0,  0,-E4,  0],
    [   0,  0,  0, E4,  0,  0],
    [   1,  0,  0,  0,  0,  0]
])

sympy.init_printing()

sympy.pprint(A.inv(method='LU'))