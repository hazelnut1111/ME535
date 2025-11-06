
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def plot_errors(h, errors):
	'''
	Plots grid spacing vs error for numeric methods on a log-log plot

	Args:
		h: 1d array of stepsizes/spacings
		errors: relative error between derivative/integral estimate and true solution
	Returns:
	'''
	plt.title("Error vs. Stepsize")
	plt.xlabel("sample spacing h")
	plt.ylabel("relative error")
	plt.loglog(h, errors)
	plt.grid()
	plt.show()





def plot_solution(t, y):
	'''
	Plots position values of Van der Pol oscillator computed using an ODE solver

	Args:
		t: 1d array of time steps
		y: 1d array of position values from Van der Pol oscillator
	Returns:
	'''

	plt.title("Numerical solution")
	plt.ylim(-2.5,2.5)
	plt.xlabel("Time t")
	plt.ylabel("Displacement y")
	plt.plot(t,y)
	plt.plot(t,y,'+')
	plt.grid()
	plt.show()


def center_diff_rad1(f,x,h):
	'''
	Compute the radius 1 centered difference estimate of the first derivative of f at x

	Args:
		f: function to evaluate
		x: value at which to estimate derivative
		h: stepsize/spacing between sample points
	Returns:
		estimate of f'(x)
	 '''
	return (f(x+h) - f(x-h))/(2*h)

def center_diff_rad2(f,x,h):
	'''
	Compute the radius 2 centered difference estimate of the first derivative of f at x

	Args:
		f: function to evaluate
		x: value at which to estimate derivative
		h: stepsize/spacing between sample points
	Returns:
		estimate of f'(x)
	'''
	
	return (-f(x+2*h)+8*f(x+h)-8*f(x-h)+f(x-2*h))/(12*h)

def compute_diff_errors(f,d_f,x,emax):
	'''
	Compute the relative error between the exact solution of f'(x) and estimates of f'(x) using a radius 2 central finite difference
	estimate with spacings 2^{-n} for n in {0,1,...,emax}.

	Args:
		f: function to evaluate
		d_f: function of the derivative of f
		x: value at which to estimate derivative
		emax: maximum exponent value
	Returns:
		h: 1d array of finite difference spacing
		errors: 1d array of relative errors between estimate of d_f(x) and exact value over various h spacings
	'''
	f=f #function
	d_f=np.cos(x) #exact derivative value
	n = emax
	e_vals1 = np.zeros(n)
	h1 = np.zeros(n)
	for i in range(0,n):
		h = 2.**(-i)
		est = center_diff_rad2(f,x,h)
		error = np.abs(d_f-est)
		e_vals1[i] = error
		h1[i]=h
		print('{:1.0e} {: 1.9f}'.format(h,e_vals1[i]))
	return(h1,e_vals1)

 
 
	
	

def writen_response_p1():
	#Based on the plot of your data, what is your estimate of the order of truncation error?
	return """ANSWER:Based on the plot, the error decreases with the second power of the sample spacing.
The estimate of the slope of the line on the log-log plot is -2, so the order of truncation error is 2.
The roundoff error becomes significant when the relative error starts to become larger at the sample spacing between 10^-4 and 10^-3. """


# def simpson38(f, xmin, xmax, p, h):
# 	'''
# 	Compute numerical quadrature value using composite Simpson's 3/8 rule

# 	Args:
# 		f: function corresponding to the integrand
# 		xmin: minimum integration range
# 		xmax: maximum integration range
# 		p: number of panels for the composition integration scheme
# 		h: panel size
# 	Returns:
# 		quadrature estimate
# 	'''
def simpson38(f, xmin, xmax, p,h):

    a=xmin
    b=xmax
    n= int((b-a)/h)
    sum= f(a)+f(b)
    for i in range(1,n):
        q = a+i*h
        if i%3==0:
            sum=sum+2*f(q)
        else:
            sum=sum+3*f(q)

    ans= (3*h/8)*sum
    return ans
    


# def compute_integ_errors(f,i_f, xmin, xmax, emax):
# 	'''
# 	Compute the relative error between the exact solution of i_f(xmax) - i_f(xmin) and estimates of i_f(xmax) - i_f(xmin) using a simpsons 3/8th
# 	quadrature with a panel count equal to 4^{n} for n in {0,1,...,emax}.

# 	Args:
# 		f: integrand function to evaluate
# 		i_f: function of the integral of f
# 		xmin: lower bound of integration range
# 		xmax: upper bound of integration range
# 		emax: maximum exponent value
# 	Returns:
# 		h: 1d array of quadrature panel sizes
# 		errors: 1d array of relative errors between estimate of i_f(xmax) - i_f(xmin) and exact value over various h spacings
# 	'''

# 	pass
def compute_integ_errors(f,i_f, xmin, xmax, emax):

    a=xmin
    b=xmax
    n=emax
    errors= [] #error save list
    h= []#h save list
    for i in range(0,n+1):
        p = 4**i
        h1=(b-a)/(3*p)
        h.append(h1)
        est =  simpson38(f, xmin, xmax, p,h1)
        exact= i_f(b)-i_f(a)
        error = np.abs(exact-est)/np.abs(exact)
        errors.append(error)
    return (h,errors)


def writen_response_p2():
	#Based on the plot of your data, what is your estimate of the order of truncation error?
	return """ANSWER:The estimate of the slope of the line on the log-log plot is -2, so the order of truncation error is 2.
The roundoff error becomes significant when the relative error starts to become larger at the sample spacing between 10^-5 and 10^-4. """


def vdp_mu(mu):
	def vdp(t, y):
		'''Van der Pol oscillator'''
		return np.array([y[1], mu*(1 - y[0]**2)*y[1] - y[0]])
	return vdp

def euler_step(f,t0,t1,y):
	"""
	Compute next value for Euler's method ODE solver
	
	Args:
		f: right-hand side function that gives rate of change of y
		y: float value of dependent variable
		t0: float initial value of independent variable
		t1: float final value of independent variable
	Returns:
		1d array of dependentant variables updated after the timestep
	"""
	h = t1 - t0
	return y + h*f(t0, y)

def rk2_step(f,t0,t1,y):
	"""
	compute next value for Euler's method ODE solver
	
	Args:
		f: right-hand side function that gives rate of change of y
		y: float value of dependent variable
		t0: float initial value of independent variable
		t1: float final value of independent variable

	Returns:
		1d array of dependentant variables updated after the timestep
	"""
	h = t1 - t0

	#compute euler estimate for half step
	y1 = y + 0.5*h*f(t0, y)
	t1 = t0 + 0.5*h
		
	#compute midstep rate estimate and take full step using midstep rate estimate 
	return y + h*f(t1, y1)

def rk4_step(f,t0,t1,y):
	"""
	compute next value for 4th-order Runge-Kutta method ODE solver
	
	Args:
		f: right-hand side function that gives rate of change of y
		y: float value of dependent variable
		t0: float initial value of independent variable
		t1: float final value of independent variable
		
	Returns:
		1d array of dependentant variables updated after the timestep
	"""
	h = t1 - t0

	#compute euler estimate for half step
	k0 = h*f(t0, y) 
	y1 = y + 0.5*k0
	t1 = t0 + 0.5*h
	k1= h*f(t1,y1)
	t2 = t0 + 0.5*h
	y2 = y + 0.5*k1
	k2 =h*f(t2,y2)
	t3 = t0 + h
	y3 = y + h*k2
	k3 = h*f(t3,y3)
	return y + (k0/6) + (k1/3) + (k2/3) +(k3/6)
	

def rk_solve(f,t,y0,order):
	"""
	Runge-Kutta solver for systems of 1st order ODEs
	
	Args:
		f: right-hand side function that gives rate of change of y
		y0: numpy array of initial float values of dependent variable
		t: numpy array of float values of independent variable
		order: int order of RK solver with allowed values [1,2,4]

	Returns:
		2D numpy array of float values of dependent variable at each time step
	"""
	if order == 1:
		method = euler_step
	elif order == 2:
		method = rk2_step
	elif order == 4:
		method = rk4_step
	else:
		raise ValueError("rk_solve `order` must be 1, 2, or 4")

	y = [y0]
	for i in range(t.size - 1):
		y.append(method(f, t[i], t[i + 1], y[i]))
	return np.array(y)

def compute_solution(f, tmin, tmax, y0, steps, order=2):
	'''
	Set up time intervals and initial conditions and compute numerical solution
	Args:
		f: function describing rates of change
		eps: float parameter in rate function
	Returns:
		t: 1D numpy float array of independent variable values
		y: 2D numpy float array of dependent variable values
	'''

	t = np.linspace(tmin, tmax, steps + 1)
	y = rk_solve(f,t,y0,order)
	return t, y

def simulate_RKF45(f, tmin, tmax, y0, steps):
	'''
	Function to set up proper arguments to perform simulation using ivp_solve.
	See scipy documentation for details.
	Note that ivp_solve may store solution data transposed from what we coded up previously.

	Return:
		t: 1D numpy array of independent variable values
		y: 2D numpy array of dependent variable values
	'''
	t_span=[tmin,tmax]
	y0=y0
	method = 'RK45'
	t_eval=np.linspace(tmin,tmax,steps)
	sol=solve_ivp(f,t_span, y0, method)
	t= sol.t
	y=	sol.y
	return t,y[1] 




def writen_response_p3():
	#3a Based on the plot of your numerical solution, describe the long-term steady-state behavior of the system.
	#3b What happens when you try this simulation with 500 steps with a Van der Pol mu of 0.1? What is the steady-state behavior with 5000 steps?
	#3c Compare and contrast the results of your RK4 and numpy's RK45 results. How many steps were used to get that agreement?
	return """ANSWER: 3(a)Because mu is small and bigger than 0, the Van Der Pol oscillator would not be a simple harmonic oscillator, it would enter a limit cycle. Based on the plot, long-term and steady-state behavior shows a stable node corresponding to a stable  2pi-periodic orbit of radius 2.  
3(b)
Based on the plot, when mu is 0.1 the Van Der Pol oscillator’s long-term and steady-state behavior has a stable node corresponding to a stable  2pi-periodic orbit of radius 2. In the plot of 500-steps,  the oscillation of the wave is between 2 and -2. When mu becomes 10, the Van Der Pol oscillator’s long-term and steady-state behavior is not stable, and its amplitude is over -+ 2, and the 5000-steps plot has more points than the 500-steps plot. 
3(c)
Based on the RKF45 plot, the oscillation of the wave is between 2 and -2 when mu is 0.1 which is similar to the RK4 plot. However, the difference between the two plots is RK4 plot has more points, the RK4 is more accurate than RKF45 in this case."""

def p1():
	'''
	1) This question involves using the center-difference formula for the first derivative from Kutz' Table 4.2 to estimate the value of df/dx for the function f(x) = sin(x) at x = 1.0. 

	For those who do not have the text near at hand, the formula (with $x$ as the independent variable and $h$ as the spacing) is:
	$ f'(x) = 1 / (12*h) * (f(x-2h) - 8*f(x-h) + 8*f(x+h) - f(x+2h))$

	Your particular tasks are as follows:

	- Compute the relative error in the derivative estimate for sample spacing values $h = 2^{-n} for n in {0,1,...,18}$.
	- Plot of the error as a function of sample spacing using the plot_errors() function.
	- Based on your plot, identify the order of truncation error and estimate the sample spacing at which the roundoff error becomes significant.
	'''

	f = np.sin
	d_f = np.cos
	x = 1.0
	emax = 18
	h,errors=compute_diff_errors(f,d_f,x,emax)
	plot_errors(h, errors)
	print(writen_response_p1())



def p2():
	'''
	2) This question involves looking at the convergence properties of a well-known method for numerical integration or quadrature known as Simpson's 3/8 rule,
	described by Eq.(4.2.6c) in the text, for the function f(x) = sin(x) between 0 and 1.

	Your particular tasks are as follows:

	- Compute the relative error in the integration estimate for panel counts of $p = 4^{n} for n in {0,1,...,9}$.
	- Plot of the error as a function of sample spacing using the plot_errors() function.
	- Based on your plot, identify the order of truncation error and estimate the sample spacing at which the roundoff error becomes significant.
	'''

	f = np.sin
	i_f = lambda x: -np.cos(x)
	xmin = 1e-12
	xmax = 1.0
	emax = 9
	h,errors=compute_integ_errors(f,i_f, xmin, xmax, emax)
	print(compute_integ_errors(f,i_f, xmin, xmax, emax))
	plot_errors(h, errors)
	print(writen_response_p2())
 

def p3():
	'''
	3) The models used by engineers often involve linear differential equations. A familiar example is the damped linear (or harmonic) oscillator: 
	$$\frac{d^2 y}{dt^2} + c \frac{dy}{dt} + y = f(t)$$ or y'' + c*y' + y = f(t)
	Linear equations are "friendly" because we can solve them analytically and they support superposition. However, the "real world" is not that
	friendly; most systems are actually nonlinear and we typically cannot write down analytic solutions. Instead, we often employ numerical methods
	to compute approximate solutions for more realistic systems modeled by nonlinear differential equations. The following problems deal with a
	classic nonlinear second-order ODE, the Van der Pol equation, which describes an oscillator with nonlinear damping:
	$$\frac{d^2 y}{dt^2} -epsilon(1-y^2) \frac{dy}{dt} + y = 0$$ or y'' - eps*(1-y**2)*y' + y = 0
	'''

	###
	# 3a) The template includes an implementation of an Runge-Kutta ODE solver `rk_solve()` that calls a function `rk2_step()` that computes a single step
	# for the second order RK method. Your initial task is to implement the function `rk4_step` to compute a single step of the 4th order Runge-Kutta
	# method. Use `rk_solve` to call your `rk4_step` function to compute a numerical solution for y(t) on the interval t=[0,100] with stepsize h = 0.1,
	# Van der Pol mu of 0.1, and initial conditions y(0)=0.25, y'(0)=0.
	###

	tmin = 0
	tmax = 100
	y0 = np.array([0.25, 0.])

	vdp = vdp_mu(0.1)
	steps = 500
	t,y = compute_solution(vdp, tmin, tmax, y0, steps)
	plot_solution(t,y)
	vdp1 = vdp_mu(10)
	steps2 = 5000
	t1,y1 = compute_solution(vdp1, tmin, tmax, y0, steps2)
	plot_solution(t1,y1)
	t2,y2=simulate_RKF45(vdp, tmin, tmax, y0, steps)
	plot_solution(t2, y2)
	t3,y3=simulate_RKF45(vdp1, tmin, tmax, y0, steps)
	plot_solution(t3, y3)
	print(writen_response_p3())

	###
	# 3b) What happens when you try this simulation with 500 steps with a Van der Pol mu of 0.1? What is the steady-state behavior with 5000 steps?
	###

	###
	# 3c) Compute a  numerical solution for the system of problem 3 using the RKF45 method  in the library function scipy.integrate.solve_ivp
	# that incorporates stepsize control. Follow the example for simulating the Lotka-Volterra system on the documentation page:
	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
	# Compare and contrast the results of the RK4 and RK45 results.
	###

if __name__ == "__main__":
	p1()
	p2()
	p3()
