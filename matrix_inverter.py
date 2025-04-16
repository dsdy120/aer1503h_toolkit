import sympy

n,e,b,a,i = sympy.symbols('n e b a i')

E1 = n*a*b*sympy.sin(i)
E2 = 0.5*n*b*sympy.cos(i)
E3 = (n*a**3*e/b)*sympy.cos(i)
E4 = 0.5*n*b
E5 = 0.5*n**2 * a
E6 = n*a**3*e/b

A = sympy.Matrix([
    [0, -E1, 0, E2, -E3, 0],
    [E1, 0, 0, 0, 0, 0],
    [0, 0, 0, E4, -E6, 0],
    [-E2, 0, -E4, 0, 0, E5],
    [E3, 0, E6, 0, 0, 0],
    [0, 0, 0, -E5, 0, 0],
])

sympy.init_printing()

sympy.pprint(A.inv(method='LU'))