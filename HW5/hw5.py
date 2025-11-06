import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import matplotlib.pyplot as plt
import scipy.special as special
import scipy.integrate as integrate
from scipy.optimize import linprog

def plot_p1(r,y):
	plt.plot(r,y.T)
	plt.title('Numerical solutions of ODE')
	plt.xlabel('r')
	plt.ylabel('y')
	plt.grid()
	plt.show()

def plot_p2(r,y,j1):
	plt.plot(r, y.T)
	plt.plot(r, j1)
	plt.plot(r, y.T - j1)
	plt.legend(['ode solution', 'special function', 'difference'])
	plt.title('Numerical solutions of ODE')
	plt.xlabel('r')
	plt.ylabel('y')
	plt.grid()
	plt.show()

def ode_rates(r, y):
	'''
	Rate functions on RHS of Problem 1 ODE as 1st-order system
	
	Args:
		r: float independent variable
		z: displacement and velocity tuple pair
	Returns:
		new displacement and velocity pair
	'''

	d, v = y
	return v, 1/(r*r)*(-r*v-(r*r-1)*d)

def library_solve(f, rmin, rmax, nr, d0, v0):
	'''
	Solve a 2nd order ode using the library function `scipy.integrate.solve_ivp`

	Args:
		f: function specifying rates for first order system
		rmin, rmax: float bounds of solution interval
		nr: int number of points at which to evaluate the solution
		d0: float initial displacement
		v0: float initial velocity
	Returns:
		2D float numpy array

	'''

	t_span = (rmin, rmax)
	t_eval = np.linspace(rmin, rmax, nr)
	sol = integrate.solve_ivp(f, t_span, (d0, v0), t_eval=t_eval)
	return sol

def root_estimates():
	'''
	Based on your numerical solution, estimate the locations 
	of the first 4 non-zero roots k1, k2, k3, k4 of your computed solution 
	(where y(k_i)==0 ).
	'''
	rmin = 1e-15
	rmax = 15
	nr = 151

	d0 = 0
	v0 = 1

	sol = library_solve(ode_rates, rmin, rmax, nr, d0, v0)
	roots = []
	for i in range(len(sol.t)-1):
		if sol.y[0][i] * sol.y[0][i+1] < 0:
			root = (sol.t[i] + sol.t[i+1]) / 2.0
			if root > 1e-10:  # exclude roots near r=0
				roots.append(root)
			if len(roots) == 4:
				break
	k1, k2, k3, k4 = roots
	return f"Root estimates: k1={k1:.1f}, k2={k2:.1f}, k3={k3:.1f}, k4={k4:.1f}"


def evaluate_bessel(rmin, rmax, nr):
	'''
	Evaluate the library function for Bessel Function J_1 on a specified grid using `scipy.integratespecial.j1`

	Args:
		rmin,rmax: float bounds of interval of interest
		nr: int number of grid points (uniformly spaced on the interval)
	Returns:
		1D float numpy array of function values sampled on the grid
	'''

	r_eval = np.linspace(rmin, rmax, nr)
	j1 = special.j1(r_eval)
	return j1

def compare_solutions():
	'''
	Compare and contrast the results from the ODE solver and the special function. 
	Do the locations of the roots coincide with the estimates from the numerical solution of the ODE?
	'''

	return """ANSWER Q2: The solutions of the ODE and the special function are very close, but not identical. 
	There are some small differences between them. The locations of the roots are the estimates 
	from the numerical solution of the ODE."""

def precise_roots(n):
	"""
	Find and evaluate the library function to determine where the zeros of a bessel function J_1 reside.

	Args:
		n: the number of desired roots for the bessel function J_1 
	Returns:
		1D float numpy array of floats specifiying the zero locations

	"""

	roots = special.jn_zeros(1, n)
	return roots

def I(ki, kj):
	'''
	Code for the integrand of problem 4 in a form suitable for scipy integration.

	Args:
		ki, kj: float  multiple of r in argument of Bessel function
	Returns:
		Integrand function f
	'''

	def integrand(r):
		'''
		Code for the integrand of problem 4 in a form suitable for scipy integration.

		Args:
			r: argument of the integrand
		Returns:
			float value of integrand
		'''

		return r * special.j1(ki*r) * special.j1(kj*r)

	return integrand

def weighted_inner(rmin, rmax, ki, kj):
	'''
	Compute the inner product (i.e. integral) of Bessel functions with weight r.

	Args:
		rmin, rmax: Float bounds of integration interval
		ki, kj: float arguments to pass to integrand
	Returns:
		q: float of the integral estimate
		error: float of the error of the integral estimate
	'''
	f = I(ki, kj)
	q, error = integrate.quad(f, rmin, rmax)
	return q, error

def inner_array(k, rmin, rmax):
	'''
	Compute the inner product of each permutation of the Bessel functions with 
	ki and kj where i,j = {1,2,3,4}

	Args:
		rmin, rmax: Float bounds of integration interval
		k: 1d float numpy array of k values to compute the inner product over
	Returns:
		2d float numpy array of inner products over the ki and kj permutations
	'''

	n=len(k)
	inner_products = np.zeros((n, n))
	for i in range(n):
		ki = k[i]
		for j in range(n):
			kj = k[j]
			q, error = weighted_inner(ki, kj, rmin, rmax)
			inner_products[i, j] = q
	return inner_products

def relation():
	'''
	Describe the array returned by `inner_array` and comment on
	the implied relationship between the functions J_1(k_i*r) and
	J_1(k_j*r).
	'''

	return """ANSWER Q3: The 'inner_array' function returns a 2D numpy array where each element (i, j) corresponds to the inner product of the Bessel functions J1(ki*r) and J1(kj*r) in the interval (rmin, rmax). 
	The relationship between these functions is that they are orthogonal to each other with respect to the weight function r, i.e., their inner product is zero unless i = j. 
	This is a general property of orthogonal functions and plays an important role in many applications."""




def random_array(n):
	'''
	Generate random symmetric array B with shape (n,n)

	Args:
		n: the length and width of the A matrix
	Returns:
		2d float numpy array of the random symmetric array B
	'''

	A = np.random.rand(n, n)
	B = 0.5 * (A + A.T)
	return B

def qr_iter(B, iters):
	'''
	Args:
		B: 2d float numpy array of the random symmetric array B
		iters: number of iterations to perform qr factorization
	Returns:
		Q: 2d float numpy array Q matrix after applying `iters` qr factorizations
		R: 2d float numpy array R matrix after applying `iters` qr factorizations
		diag: 1d float numpy array of the diagonal elements of Q @ R
	'''

	for i in range(iters):
		Q, R = np.linalg.qr(B)
		B = np.dot(R, Q)
	diag = np.diag(B)
	return Q, R, diag
		
def purpose():
	'''
	What is the outcome of the iterated QR-factorization scheme?
	'''

	return """ANSWER Q4: The iterated QR-factorization scheme converges to the diagonal of the original matrix. The 
	non-zero entries in the R matrix indicate the off-diagonal elements that remain after each iteration."""

def my_linprog(equality=True):
	'''
	Refer to PDF for problem description
	'''

	obj = [-1, -2]
	lhs_ineq = [[1, 1], [-4, 5], [1, -2]]
	rhs_ineq = [10, 10, 2]
	lhs_eq = [[-1, 5]]
	rhs_eq = [15]
	bnd = [(0, None), (0, None)]
	
	if equality:
		opt = linprog(obj, A_eq=lhs_eq, b_eq=rhs_eq, A_ub=lhs_ineq, b_ub=rhs_ineq, bounds=bnd)
	else:
		opt = linprog(obj, A_ub=lhs_ineq, b_ub=rhs_ineq, bounds=bnd)
	
	return opt

def p1():
	rmin = 1e-15
	rmax = 15
	nr = 151

	d0 = 0
	v0 = 1

	sol = library_solve(ode_rates, rmin, rmax, nr, d0, v0)

	plot_p1(sol.t, sol.y[0])

	print(root_estimates())

def p2():
	rmin = 1e-15
	rmax = 15
	nr = 151

	d0 = 0
	v0 = 1

	sol = library_solve(ode_rates, rmin, rmax, nr, d0, v0)
	j1 = evaluate_bessel(rmin, rmax, nr)

	plot_p2(sol.t, sol.y[0], j1)

	print(compare_solutions())

def p3():
	k = precise_roots(4)
	rmin = 0
	rmax = 1
	arr = inner_array(k, rmin, rmax)

	print(relation())


def p4():
	n = 4
	iters = 10

	B = random_array(n)
	Q, R, diag = qr_iter(B, iters)

	plt.spy(R)
	plt.show()

	eigs = np.linalg.eigvals(B)
	diff = np.sort(eigs) - np.sort(diag)
	plt.plot(np.sort(eigs),'o', np.sort(diag), 'x', diff)
	plt.legend(["Eigenvalues", "Diagonal", "Difference"])
	plt.show()

	print(purpose())

def p5():
	
	# your solutions should match the following
	# 	With equality contraint:
	# 	Optimum  -14.17 at (x,y) =  [5.83 4.17]
	# 	Without equality contraint:
	# 	Optimum  -15.56 at (x,y) =  [4.44 5.56]

	np.set_printoptions(precision=2)
	opt = my_linprog(equality=True)
	print("With equality contraint:")
	print("Optimum ", opt.fun, "at (x,y) = ", opt.x)

	opt = my_linprog(equality=False)
	print("Without equality contraint:")
	print("Optimum ", opt.fun, "at (x,y) = ", opt.x)

if __name__ == '__main__':
	p1()
	p2()
	p3()
	p4()
	p5()