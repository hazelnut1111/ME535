import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import matplotlib.pyplot as plt

def plot_p1(v):
	'''
	Plots 2d v array for problem 1
	'''
	m, n = v.shape
	#note: swap positions of elements in v so that the triangle
	#is in the bottom left and +x/+y are in the right/up direction
	v = np.flipud(np.fliplr(v))
	x = np.linspace(0,1,m)
	y = np.linspace(0,1,n)
	plt.contourf(x, y, v)
	plt.show()

def plot_p2(x, t, u):
	'''
	Plots 2d u array for problem 2
	'''

	X,T = np.meshgrid(x,t)
	fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
	surf = ax.plot_surface(X, T, u, cmap=cm.coolwarm, linewidth=0, antialiased=False)

	ax.set_zlim(-1.01, 1.01)
	ax.zaxis.set_major_formatter('{x:.02f}')
	fig.colorbar(surf, shrink=0.5, aspect=5)
	plt.show()

def initialize_1(u, h):
	'''
	Initialize an array to account for BCs and geometry for problem 1.

	Args:
		u: 2D numpy array corresponding to values sampled on a regular 2D grid
		h: float grid-spacing
	Returns:
		2D numpy array corresponding to the input with boundary values on the boundary/exterior or the triangular domain.
	'''

	pass

def update(u):
	'''
	Update values within the triangular domain based on the 5-point averaging stencil

	Args:
		u: 2D numpy array of input values (with appropriate boundary conditions outside the triangular domain)
	Returns:
		2D numpy array of with values updated within the triangular domain
	'''

	pass

def iter_solve(u, tol, max_iters):
	'''
	Solve the Laplace equation by iterating with 5-point averaging stencil (by calling update() )   
	
	Args:
		u: 2D numpy array to update
		tol: float maximal allowable change in an updated value
		max_iters: int limit on number of update iterations

	Return:
		2D numpy array with final solution
		number of iterations to achieve specified tolerance
	'''

	pass

def displacement(x):
	'''
	Initial displacement function for problem 2

	Args:
		x: 1d numpy float array of positions on uniform 1d grid
	Returns:
		1d numpy float array of displacements
	'''

	return np.maximum(0, 1 - np.abs(x/10.0))

def initialize_2(f,x,u):
	'''
	Set up the solution array and establish initial conditions

	Args:
		f: function that defines initial displacement
		x: 1d numpy float array of positions on uniform 1d grid
		u: nt x nx float numpy array
	Returns:
		nt x nx float numpy array with initial displacement stored for first two time steps
	'''

	pass

def single_step(x0,x1,dx,dt):
	'''
	Compute a single update step for central difference 
	discretization of wave equation.

	Args:
		x0: 1D numpy float array of previous displacements
		x1: 1D numpy float array of current displacements
		dx: float grid-spacing
		dt: float timestep
	Returns:
		x2: 1d array of updated displacements
	'''

	pass

def step_solve(x,t,u,dx,dt):
	'''
	Compute the solution array for problem 2

	Args:
		x: 1d numpy float array of positions on uniform 1d grid
		t: 1d numpy float array of time steps on uniform 1d grid
		u: nt x nx float numpy array
		dx: float grid-spacing
		dt: float timestep
	Retuns:
		2D float numpy array of computed solution values

	'''

	pass

def p1():
	#refer to PDF for problem description

	n = 10
	h = 1/(n-1)

	tol = 0.01
	max_iters = 30

	u = np.zeros([n,n])
	# u = initialize_1(u, h)
	# v, iters = iter_solve(u, tol, max_iters)
	# plot_p1(v)

def p2():
	#refer to PDF for problem description

	xmin = -100
	xmax = 100
	nx = 501
	dx = (xmax - xmin)/(nx - 1)

	tmin = 0.0
	tmax = 35.0

	for s in [0.999,1.001]:
		dt = s*dx
		nt = int((tmax - tmin)//dt) + 1

		x = np.linspace(xmin,xmax,nx)
		t = np.linspace(tmin,tmax,nt)
		u = np.zeros((nt,nx))

		# u = initialize_2(displacement,x,u)
		# u = step_solve(x,t,u,dx,dt)
		# plot_p2(x,t,u)

if __name__ == '__main__':
	p1()
	p2()
