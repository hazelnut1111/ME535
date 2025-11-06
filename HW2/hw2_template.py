

import numpy as np
from numpy.testing import assert_allclose

np.set_printoptions(precision = 4)

def row_op(A,c,i,j):

    """
	perform elementary row operations
	if i == j, multiply row i by constant c
	if i != j, add c*(row i) to row j
	
	Args:
		A: 2D numpy array representing a matrix
		c: multiplicative constant
		i,j : row indices
	"""

    if i == j:
        A[i]=A[i]*c 
    else:
        A[j]=A[i]*c+A[j]
    return A




def swap_rows(A,i,j):
    """
	perform elementary row operation to swap rows i and j
	
	Args:
		A: 2D numpy array representing a matrix
		i,j: integer row indices
	"""
    a = np.copy(A[i])
    A[i] = A[j]
    A[j] = a
    return A


	

def dominant_eigen_iteration(A, u0, tol, max_iters):
    """
	compute dominant eigenvector and eigenvalue for square matrix
	
	Args:
		A:  nxn numpy array representing a matrix
		u0: initial estimate of eigenvector (1D numpy array)
		tol: float relative error termination criterion
		max_iters: integer iteration count termination criterion
	Returns:
		lambda: float dominant eigenvalue
		v: dominant eigenvector 1d float numpy array 
	"""
    n=0
    save_list=[0]
    for i in range(max_iters):
        # normalizing u0
        v = u0/np.linalg.norm(u0)
        # calculate non-normalizing eigenvector
        u0 = A.dot(v)
        # calculate eigenvalue lambda
        lam = v.dot(u0)
        #save lambda in the list
        save_list.append(lam)
        #iteration count
        n=n+1
        #stop iterating
        if (np.abs(save_list[i]-save_list[i-1])<tol):
            break
        # normalizing final eigenvector
    v = u0 / np.linalg.norm(u0)
    print('max eigenvalue:',lam)
    print('eigenvector:',v)
    print('iteration:',n)
        
    return v, lam 
# A = np.array([[1, 4,5,1],[2,2,2,2],[3,1,2,5],[6,1,2,7]])
# u0 = np.array([1,2,5,1])
# max_iters = 10000
# tol=1e-10
# print(dominant_eigen_iteration(A, u0, tol, max_iters))

	

def recessive_eigen_iteration(A, u0, tol, max_iters):
    """
	compute recessive eigenvector and eigenvalue for square matrix
	
	Args:
		A:  nxn numpy array representing a matrix
		u0: initial estimate of eigenvector (1D numpy array)
		tol: float relative error termination criterion
		max_iters: integer iteration count termination criterion
	Returns:
		lambda: float recessive eigenvalue
		v: recessive eigenvector 1d float numpy array 
	"""
    n=0
    save_list=[0]
    for i in range(max_iters):
        # normalizing u0
        v = u0/np.linalg.norm(u0)
        # Solve A v^{(k)} = u^{(k−1)}
        u0 = np.linalg.solve(A, v)
        # calculate eigenvalue lambda
        lam = 1/v.dot(u0)
        #save lambda in the list
        save_list.append(lam)
        #iteration count
        n=n+1
        #stop iterating
        if (np.abs(save_list[i]-save_list[i-1])<tol):
            break
        # normalizing final eigenvector
    v = u0 / np.linalg.norm(u0)
    print('min eigenvalue:',lam)
    print('eigenvector:',v)
    print('iteration:',n)
        
    return v, lam  
# print(recessive_eigen_iteration(A, u0, tol, max_iters))


	

def condition(A,u0, tol, max_iters):
    '''
	Compute numerical estimate of condition number of a matrix based on eigenspectrum
	Args:
		A: 2D numpy array representing the matrix
		u0: 1D numpy array that serves as initial guess for eignvector iteration
		tol: float residual for termination
		max_iters: int bound on number of iterations
	Returns:
		estimate of condition number
	'''
    v1,lam1 = dominant_eigen_iteration(A, u0, tol, max_iters)
    print(lam1)
    v2,lam2 = recessive_eigen_iteration(A, u0, tol, max_iters)
    print(lam2)
    estimate = np.abs(lam1-lam2)
    return estimate



def component_of_along(v, u):
    '''
	Compute component of vector v along direction of vector u

	Args:
		v,u: 1d numpy arrays
	Returns:
		1d numpy array representing the component of v along u
	'''


    v = np.array(v)
    u = np.array(u)
  # e = the unit vector of u
  # e = u /np.linalg.norm(u)
  # the component of v along u = |v|*cos(theta) multiplied by unit vector of u
    m =  (np.dot(v,u)/pow(np.linalg.norm(u),2))*u
    return m


	

def reflect(v,u):
    '''
	Compute reflection of vector v across mirror hyperplane with normal vector u

	Args:
		v,u: 1d numpy arrays
	Returns:
		1d numpy array representing the refelction of v
	'''
    v = np.array(v)
    u = np.array(u)
    # ut is the transpose of u
    ut= np.transpose(u)
    # unit vector of u and ut
    u1= u /np.linalg.norm(u)
    ut1=ut /np.linalg.norm(ut)
    # the reflection of vector v is equal to (v-2*u1*ut1*v)
    v1 = v-2*u1*ut1*v
    return v1



def reflect_to_e0(u):
    '''
	Compute the matrix that rotates a given vector to the e0 direction

	Args:
		u: 1D numpy array representing vector to rotate
	Returns:
		reflection: 2D numpy array representing the rotation matrix
	'''

    size = len(u)
    index = 0
    e0 = np.array([1.0 if i == index else 0.0 for i in range(size)])
    v1 = np.array(e0)
    v2 = np.array(u)

         #normalizing
    n1 = v1 / np.linalg.norm(v1)
    n2 = v2 / np.linalg.norm(v2)
    
    dot_product=np.dot(n1,n2)
        # rotation by angle
    angle = np.arccos(dot_product)

    I = np.identity(size)
        # I find rotation matrix in n is
        # R=I+(𝑛2𝑛𝑇1−𝑛1𝑛𝑇2)sin(𝑎)+(𝑛1𝑛𝑇1+𝑛2𝑛𝑇2)(cos(𝑎)−1)
    R = I + ( np.outer(n2,n1) - np.outer(n1,n2) ) * np.sin(angle) + ( np.outer(n1,n1) + np.outer(n2,n2) ) * (np.cos(angle)-1)
    return R
# u = [1,2,3,4,5,6,7]
# print(reflect_to_e0(u))








def Householder(A):
    '''
	Compute QR0 partial matrix factorization based on Householder reflection

	Args:
		A: 2D numpy array representing matrix to be factored
	Returns
		Q: 2D float numpy array representing othrogonal factor
		R0: 2D float numpy array representing whose first column is e0
	'''
    m = A.shape[0]
    Q = np.identity(m)
    R0 = np.copy(A)
    for i in range(1):
        x = R0[i:, i]
        e = np.zeros_like(x)
        e[0] = np.linalg.norm(x)
        u = x - e
        v = u / np.linalg.norm(u)
        H = np.identity(m)
        H[i:, i:] =H[i:, i:]- 2.0 * np.outer(v, v)
        # R=H(n-1)*...*H(1)*A
        R0 = H.dot(R0)  
        ## Q=H(n-1)*...*H(1)
        Q = Q.dot(H)  
    return Q, R0




# A = np.array( [   [1.0, 2.0,1.0],
#                     [2.0, 1.0,2.0],
#                       ])

# q, r = Householder(A)
# print("Householder")
# print("Q: \n", q)
# print("R: \n", r)
# print(q.dot(r))


# A = np.array( [   [1.0, 2.0],
#                     [2.0, 1.0],
#                       ])

# q, r = Householder(A)
# print("Householder")
# print("Q: \n", q)
# print("R: \n", r)



def QR_Householder(A):
    '''
	Compute QR matrix factorization based on Householder reflection

	Args:
		A: 2D numpy array representing matrix to be factored
	Returns
		Q: 2D float numpy array representing othrogonal factor
		R: 2D float numpy array representing upper triangular factor
	'''
    m = A.shape[0]
    Q = np.identity(m)
    R = np.copy(A)
    for i in range(m-1):
        x = R[i:, i]
        e = np.zeros_like(x)
        e[0] = np.linalg.norm(x)
        u = x - e
        v = u / np.linalg.norm(u)
        H = np.identity(m)
        H[i:, i:] =H[i:, i:]- 2.0 * np.outer(v, v)
        # R=H(n-1)*...*H(1)*A
        R = H.dot(R)  
        ## Q=H(n-1)*...*H(1)
        Q = Q.dot(H)  
    return Q, R

# A = np.array( [   [0.5, 1.0],
#                     [2.0, 1.0],
#                     [2.0, 1.0]
#                       ])

# q, r = QR_Householder(A)
# print("Householder")
# print("Q: \n", q)
# print("R: \n", r)
# print(q.dot(r))



def p1():
	'''
	1) This problem involves writing and using an alternate implementation of elemetary row operations.
	The implementation in the Ch. 2 notebook uses a `for` loop to individually update each entry in
	the row being altered. Here you will write an alternate version supported by `numpy` that provides
	a way to refer to a portion of an array including multiple elements including an entire (or partial)
	row or column.

	The related keywords in `numpy include "views", "slice", and "copy". Please read the materials at the
	following link that provides a nice description of the details:
	https://jakevdp.github.io/PythonDataScienceHandbook/02.02-the-basics-of-numpy-arrays.html
	'''

	'''
	1a) Write a new version of `row_op` that operates row-wise instead of element-wise.
	'''

	mat = np.eye(3)
	row_op(mat,2,1,1)
	ans = np.array([[1,0,0],[0,2,0],[0,0,1]], dtype = np.float64)
	assert_allclose(mat, ans, err_msg='failed first row op test')

	mat = np.eye(3)
	row_op(mat,3,1,2)
	ans = np.array([[1,0,0],[0,1,0],[0,3,1]], dtype = np.float64)
	assert_allclose(mat, ans, err_msg='failed second row op test')

	'''
	1b) Previously we skipped one of the elementary row operations; i.e. swapping rows. Now is the time to fill in that gap. Write code for a row-wise implementation of `swap_rows(A, i0, i1)` that exchanges rows `i0` and `i1` in the array `A`.
	'''

	mat = np.eye(3)
	swap_rows(mat,0,1)
	ans = np.array([[0,1,0],[1,0,0],[0,0,1]], dtype = np.float64)
	assert_allclose(mat, ans, err_msg='failed row swap test')

def p2():
	'''
	2a) Implement the matrix iteation scheme for computing the "dominant" eigenvalue/eigenvector (associated with the eigenvalue with largest magnitude) of a square matrix.
	Here is one verion of a simple pseudo-code (that needs termination conditions)
	Choose an initial vector u0 such that ∥u0∥ = 1
	for k = 1, 2, . . . do
		v^{(k)} = A u^{(k−1)}
		u^{(k)} = v^{(k)}/∥v^{(k)}∥
	end
	'''

	'''
	2b) Implement the "inverse matrix iteration" scheme for computing the "recessive" eigenvalue/eigenvector pair (associated with the eigenvalue with smallest magnitude) of a square matrix. 
	Here is a simple pseudo-code (again missing termination conditions)

	Choose an initial vector u0 such that ∥u0∥ = 1
	for k = 1, 2, . . . do
		Solve A v^{(k)} = u^{(k−1)}
		u^{(k)} = v^{(k)}/∥v^{(k)}∥
	end

	Use `numply.linalg.solve` in your implementation.
	'''

	'''
	2c) Use the functions you implemented in parts a and b to create a function that computes an estimate of the condition number of a square matrix.
	'''
	pass

def p3():
	'''
	3) This problem involves functions to support QR factorization.
	'''

	'''
	3a) Implement a function to compute the component of vector v along the direction of vector u.
	'''

	'''
	3b) In class, we discussed QR factorization based on Gram Scmidt orthogonalization. In that approach, an orthogonal basis is constructed by subtracting from each new candidate basis vector the components along the direction of each already-computed entry in the (numerically) orthogonal set. That approach can run into precision issues (because candidate basis vectors that lie close to the space spanned by the existing basis vectors can lead to catastrophic cancellation when components are subtracted). That issue led to the creation of other approaches. The one that is the focus of this problem involves Householder reflections. The basic idea is to compute the vector that arises when an input vector v is reflected about the plane normal to a specified vector u.

	Fill in code below to implement a function to compute the Householder reflection. This should be pretty straightforward if you usethe function component_of_along that you implemented in 3a.
	'''

	'''
	3c) Householder made good use of this basic reflection operation. First, he noted that the reflection operation corresponds to multiplication of the input vector v by an orthogonal matrix Q_u. Think about why that is true! Then, he figured out how to pick the mirror normal that would produce a reflected vector e0 that lies along the first coordinate axis; i.e. Q_u(v) = norm(v)*e0. Convince yourself that the choice u = c*(v - norm(v)*e0) works for any choice of the constant c. In particular, choose c=1 so u = v - norm(v)*e0

	NOTE: The slightly more complicated choice u = v - sign(v[0])*norm(v)*e0 prevents the possibility of catastrophic cancellation errors.

	Insert code below to produce the orthogonal matrix that reflects a given input vector to become a multiple of e0.
	'''

	'''
	3d) Use the function you implemented for part 3c to implement a function that does the first step of QR factorization.
	'''

	pass

def p4():
	'''
	4) OPTIONAL: Use the functions you wrote for problem 3 to implement QR factorization 
	based on Householder reflections.
	'''
	pass

if __name__ == '__main__':
	p1()
	p2()
	p3()
	p4()
