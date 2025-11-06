---
theme: classic #gaia
paginate: true
# class:
#  - lead
size: 16:9

marp: true
footer: "ME 535 Winter 2023, Prof. Storti, UW-Seattle" 

---


# Linear Systems (a.k.a. systems of linear equations)

to accompany Ch. 2 of __Data-Driven Modeling & Scientific Computation__ by J. Nathan Kutz

<!-- 
Some notes here that might be useful (but only in Powerpoint or HTML prsentation).
-->

---

- __linear__ $\iff$ each term involves an unknown with a constant coefficient
- __affine__ $\iff$ includes constant terms
- Start with $n$ equations in $n$ unknowns $x_0, x_1, \ldots , x_{n-1}$:

$$\begin{aligned}
a_{0,0} x_0 &+ a_{0,1} x_1 &+ \ldots &+ a_{0,n-1} x_{n-1} &= b_0& \\
a_{1,0} x_0 &+ a_{1,1} x_1 &+ \ldots &+ a_{1,n-1} x_{n-1} &= b_1& \\
\vdots & &\ddots &  &\vdots &\\
a_{n-2,0} x_0 &+ a_{n-2,1} x_1 &+ \ldots &+ a_{n-2,n-1} x_{n-1} &= b_{n-2}& \\
a_{n-1,0} x_0 &+ a_{n-1,1} x_1 &+ \ldots &+ a_{n-1,n-1} x_{n-1} &= b_{n-1}&
\end{aligned}$$
- Left-hand side (LHS) is a __linear system__; all together, a set of __affine__  equations

---

Nobody wants to write (or typeset) that over and over so...
- Collect variables and RHS as __vectors__ (1D arrays) 	
	- Typical elements $x_j, b_i$
- Collect coefficients as __matrix__ (2D array)
  -  Typical element $a_{i,j} \; :  \quad i \iff$ row; $j \iff$ column

-Abbreviations: 
- $\sum_{j=0}^{n-1} a_{i,j} x_j = b_i$
  - Einstein Summation Convention $\implies a_{i,j} x_j = b_i$ (sum on repeated index)
  - Coordinate free notation: $A x = b$

---

## Section 2.1 - Direct Solution Methods for $A x=b$

Start with the original solution method: __Gaussian elimination__. 

>Note likely historical misattribution. According to Wikipedia, first introduced in "Western" math by Newton, but the first used in China possibly 2 millenia ago.

- Basic matrix operations provide opportunites to use `for()` loops
- Take advantage of opportunities to develop our execution control skills.

- __Basic idea__: 
  - Create __augmented matrix__ $[A | b ]$ (attach RHS as extra column) 
  - Perform elementary row operations until the matrix becomes the identity. 
  - At the end of that process, augmented column contains the solution $x$. 
  - Conceptually simple, but not most efficient...

---

Approach: 
- Code up simpler things and work our way up to a full solver. Start by:
-  Importing `numpy` (to get access to `array` capabilities)
-  Creating a sample matrix (2D array) and vector (array with a single row or column).

Example in text involves matrix $a=[[1,1,1],[1,2,4],[1,3,9]]$. 
There are many ways to create a corresponding array in python. The "manual" approach (where you enter the desired values from the keyboard) is straightforward, but perhaps not terribly instructive. 

---

Manually construct array and print array followed by its transpose.
```python
import numpy as np #import numpy for access to array capabilities
a = np.array([[1,1,1],[1,2,4],[1,3,9]])
print('a=\n',a, "\ntranspose of a=\n", a.T)
```
a=
 [[1 1 1]
 [1 2 4]
 [1 3 9]] 
transpose of a=
 [[1 1 1]
 [1 2 3]
 [1 4 9]]

 ---

 Transpose `a.T` suggests algorithmic approach (generalizes to different array sizes). 
Each row in the transpose contains powers of a sequence of integers starting with 1. 
Code below sets the matrix size and creates a corresponding array of zeros. 
```python
n = 3 #set the size of the square matrix
a_transpose = np.zeros([n,n]) # create a 3x3 array initialized with zeros
# reassign the entries according to a_transpose[i,j] = 1+i+j**2
for i in range(n):
    for j in range(n):
        a_transpose[i,j] = (j+1)**i
a = a_transpose.T
print(a)
```
[[1. 1. 1.]
 [1. 2. 4.]
 [1. 3. 9.]]

---

Define function to construct $n \times n$ matrix of this type:

```python
def make_a(n):
    """
    construct n x n matrix with columns containing powers of natural numbers
	Inputs: n, integer number of rows/columns
	Returns: 2D array
    """
    a = np.zeros([n,n]) #create array of size needed
    for i in range(n): #loop over rows
        for j in range(n): #loop over columns
            a[i,j] = (j+1)**i #assign value
    return a.transpose()

a = make_a(3); print('a=\n',a) #call function and print result
```
a=
 [[1. 1. 1.]
 [1. 2. 4.]
 [1. 3. 9.]]

 ---

 Construct array corresponding to RHS: $b = [1,-1,1]$:
 ```python
b = np.array([1,-1,1])
print('b=',b)
 ```
b= [ 1 -1  1]
Again construct function to produce arrays of given size:
```python
def make_b(n):
    a = np.zeros(n)
    for i in range(n):
        a[i]=(-1)**i
    return a

print('b(',n,')=', make_b(7))
```
b(7)= [ 1. -1.  1. -1.  1. -1.  1.]

---

> Aside: $b$ looks like a row-vector. Column version can be produced in various ways including:
```python
print('transpose version:\n', np.transpose([b]))
print('newaxis version:\n',b[:, np.newaxis])
```
Test them out!

---

## Starting toward a solver
### Construct augmented matrix by appending b as an extra column of a
```python
a = make_a(3)
b = make_b(3)
auga = np.column_stack([a,b])
# alternate using np.c_ (check docs for details)
# auga= np.c_[a,b]
print('augmented matrix:\n',auga)
```
augmented matrix:
 [[ 1.  1.  1.  1.]
 [ 1.  2.  4. -1.]
 [ 1.  3.  9.  1.]]

 ---

 Create function for specified size:
 ```python
def augment(A,b):
    '''
    Construct augmented matrix
    Inputs:
        A: 2-dim numpy array of shape (n,n)
        b: 1- numpy array of length n
    Returns:
        n x (n+1) numpy array
    '''
    m,n = A.shape
    if m!=n:
        print("Input A is not square.")
    return np.column_stack[A,b]
auga = augment(a,b)
print('The augmented matrix is:\n', auga)
 ```
 The augmented matrix is: [[ 1.  1.  1.  1.]
 [ 1.  2.  4. -1.]
 [ 1.  3.  9.  1.]]

 ---

 ## Simplify using __elementary row operations__

 - Swap rows
 - Multiply row by a non-zero constant
 - Add multiple of one row to another

__Elementary row ops leave solution unchanged__
- After implementation, apply elementary row ops to systematically transform/simplify until we can readily solve
- Be a little clever and implement them all in a single function

---

```python
def row_op(A,c,i1,i2):
    """
    perform elementary row operations on 2D array A
    if i1==i2, multiply row by constant c
    if i1!=i2, add c*(row i1) to row i2
    Args:
        A: 2D numpy array representing a matrix
        c: float multiplicative constant
        i1,i2: int row indices
    """
    m,n = A.shape #number of rows and columns
    if  i1<0 or i2<0 or i1>=m or i2>=m:
        print("WARNING: Invalid index specifications. Each index value i must satisfy 0<=i<#rows.")
    if i1==i2: #repeated index -> multiply row by constant
        for j in range(n):
            A[i1,j] *= c #Equiv. to A[i1,j] = c*A[i1,j]
    else: # add c*row i1 to row i2
        for j in range(n):
            A[i2,j] += c * A[i1,j] # Equiv. to A[i2,j] = A[i2,j] + c*A[i1,j]
    return
```

---

## Simplification scheme
Use row ops to obtain __triangular__ form (entries below main diagonal are 0). How?
- Choose a diagonal element `A[p,p]` as a "pivot".
- Divide row `p` by pivot (to produce `A[p,p]` = 1 at pivot on diagonal)
- Perform row operations to zero out each coefficient beneath it in column `p`
    - Subtract `A[i,p]` times pivot row `p` from row `i`.
- If preceding entries in a row are already zero, they remain zero during  elementary row ops.
- Proceed systematically across the columns:
    - Zero out below diagonal to produce upper triangular system.

---
__Implement__: Write `cancel_below_diagonal()` to zero out elements in lower triangle
```python
def cancel_below_diagonal(A, pivot):
    """
    insert docstring here 
    """
    SMALL_VALUE = 1E-8 # check for possible overflow
    m,n = A.shape
    if  pivot<0 or pivot>=m:
        print("WARNING: Invalid index specification. Index value pivot must satisfy 0<=pivot<#rows.")
    if abs(A[pivot,pivot]) < SMALL_VALUE:
        print("WARNING: Division by near-zero pivot value.")
    else:
        # row_op(A,1./A[pivot,pivot], pivot, pivot) #divide by pivot value so value at pivot position becomes 1
        for i in range(pivot+1,m):
            row_op(A,-A[i,pivot]/A[pivot,pivot],pivot,i)   
```

---

```python
#test the function by zeroing the first column below the diagonal
auga = augment(a,b)
cancel_below_diagonal(auga,0)
print("Original matrix A = \n", a)
print("Augmented matrix after processing 1st pivot = \n", auga)
```
Original matrix A = 
 [[1. 1. 1.]
 [1. 2. 4.]
 [1. 3. 9.]]
Augmented matrix after processing 1st pivot = 
 [[ 1.  1.  1.  1.]
 [ 0.  1.  3. -2.]
 [ 0.  2.  8.  0.]]

 ---

 ```python
#continue by zeroing the second column below the diagonal
cancel_below_diagonal(auga,1)
print("Augmented matrix after processing 2nd pivot = \n", auga)
 ```
 Augmented matrix after processing 2nd pivot = 
 [[ 1.  1.  1.  1.]
 [ 0.  1.  3. -2.]
 [ 0.  0.  2.  4.]]

 Successfully achieve upper triangular form (with 0 below main diagonal `A[i,i]`)

 ---

 Define function `upper_tri()` to triangularize in single call:
 ```python
def upper_tri(A):
    m,n = A.shape
    for i in range(m):
        cancel_below_diagonal(A,i)

auga = augment(a,b)
upper_tri(auga)
print("Inspect to verify that augmented array has been triangularized:\n", auga)
 ```
Inspect to verify that augmented array has been triangularized:
 [[ 1.  1.  1.  1.]
 [ 0.  1.  3. -2.]
 [ 0.  0.  2.  4.]]

 ---
## Backsolve/backsubstitution
- Given triangular matrix, the last row corresponds to a linear equation in 1 variable that can be solved immediately: $x_2 = A_{3,4}/A_{3,3} = 4./2. = 2.$
- Back-substitute value of $x_2$ into rows above.
- Next row up gives equation to solve for $x_1$.
- Plug in above; solve for $x_0$ to complete solution. 

Are you done at that point?
Check your result. 
How? Compute the residual $r = Ax-b$; should be close to $\vec{0}$.

---

Implement a `back_sub()` function to execute the back substitution process:

```python
def back_sub(augU):
    """
    Insert suitable docstring here.
    """
    m,n = augU.shape
    x = np.zeros(m)
    for i in range(m):
        x[m-1-i]= augU[m-1-i,-1] #Initialize solution entry with value from RHS of tri. system
        for j in range(m-i,m): #For each entry of the row right of the main diagonal
            x[m-i-1] -= augU[m-i-1,j]*x[j] #Subtract coeff. * (known/larger-index entry in solution)
        x[m-1-i] /= augU[m-1-i, m-1-i] #Divide by pivot to get the new entry in the solution
    return x  
```

---
Test functions for triangularizing and back-substituting:

```python
a = make_a(3)
b = make_b(3)
auga = augment(a,b)
upper_tri(auga)
print("Triangular aumented matrix:\n", auga)
soln = back_sub(auga)
print("Solution obtained by back-substitution:", soln)
print("Check that solution satisfies the equations:\nResidual = ", np.dot(a, soln) - b)
```
Triangular aumented matrix:
 [[ 1.  1.  1.  1.]
 [ 0.  1.  3. -2.]
 [ 0.  0.  2.  4.]]
Solution obtained by back-substitution: [ 7. -8.  2.]
Check that solution satisfies the equations:
Residual =  [0. 0. 0.]

---

Modify test problem with more pivots not equal to 1:
```python
a = np.array([[3., 3., 3.], [2., 4., 8.], [1., 3., 9.]])
b = np.array([3., -2., 1.])
auga = augment(a,b)
upper_tri(auga)
print("Triangular aumented matrix:\n", auga)
soln = back_sub(auga)
print("Solution obtained by back-substitution:", soln)
residual = np.dot(a, soln) - b
print("Residual = ", residual, ", Norm of the residual = ", np.linalg.norm(residual))
```
Triangular aumented matrix:
 [[ 3.  3.  3.  3.]
 [ 0.  2.  6. -4.]
 [ 0.  0.  2.  4.]]
Solution obtained by back-substitution: [ 7. -8.  2.]
Residual =  [0. 0. 0.], Norm of the residual =  0.0

---

```python
# check that Ax agrees with b to within a threshold
# print the entries in Ax-b
print("Ax-b = ", a.dot(soln)-b)
# use numpy's `allclose` function to check for numerical agreement
# to within a specified absolute and/or relative tolerance
print("Solution checks? : ", np.allclose(np.dot(a,soln),b, atol = 1e-10))
```
Ax-b =  [0. 0. 0.]
Solution checks? :  True

---

This approach basically works (at least in the cases we have seen so far), but the whole process needs to run again to solve with a different right-hand side `b`. Instead of repeating the computation, it is more convenient to rephrase the problem in terms of matrix factorization. The particular factorization of interest here is called the __LU Factorization__ because it involves rewriting $A$ as the product of a lower-triangular matrix $L$ and the upper-triangular $U$ that was computed above by row reduction.

Once the factorization is computed, the solution for any right-hand side can be found by solving the 2 triangular systems $L y = b$ and then $U x = y$.

>Note that the 2 equations together are equivalent to the original system:
<br>$Ly=b \land Ux=y \implies L(Ux)=b \iff (LU)x=b \iff Ax=b$.

__Next goal__: Write function to compute the $LU$ factorization and to solve given $L,U,b$. 

---

First, some context:

Previously discussed solution by __Gaussian Elimination__:

- Perform elementary row ops to triangularize
- Solve triangular system by back substitution
- Bascally works (at least for simple cases seen so far)
- What if we need to solve again with different $b$ on RHS?

Instead of starting over, "record" row ops for triangularization.

Change our perspective from row ops to matrix multiplication:
- What matrix multiplication is needed to undo the row operation?
- Keep both the updated matrix and the matrix that undoes the update (to preserve the initial matrix as the product)

---

Back to example:
```python
a = make_a(3)
b = make_b(3)
auga = augment(a,b)
print(auga)
```
 [[ 1.  1.  1.  1.]
 [ 1.  2.  4. -1.]
 [ 1.  3.  9.  1.]]

 ---

First triangularization step: Subtract Row 0 from Rows 1 and 2 to zero out the subdiagonal entries in Column 0.
<br>What operation undoes the cancellation? Add Row 0 to Row 1 and Row 2.
<br>How do we write that in the language of linear algebra?

When multiply matrix $A$ by column vector $x$, each entry in the output is obtained by multiplying a row by the entries in  $x$ and summing.

Net result: Multiply a matrix by a column vector (on the right) to produce a linear combinations of the columns of  $A$ with coefficients in $x$.

---

Here we want linear combinations of rows, so we transpose: 

- Pre-multiply the matrix by a row of coefficients to produce the corresponding linear combination of the rows.
- Stack the rows to form a matrix.

$$\begin{aligned}
\begin{Bmatrix}
R_1 &\leftarrow & R_1 \\
R_2 &\leftarrow & 1*R_1 + R_2 \\
R_3 &\leftarrow & 1*R_1 + R_3
\end{Bmatrix} \iff
\begin{pmatrix}
1 & 0 & 0 \\
1 & 1 & 0 \\
1 & 0 & 1
\end{pmatrix}
A
\end{aligned}$$

---

Do the factorization described to clear the first column below the diagonal:

$$\begin{aligned}
A = 
\begin{pmatrix}
    1 & 1 & 1 \\
    1 & 2 & 4 \\
    1 & 3 & 9
\end{pmatrix}
= 
\begin{pmatrix}
    1 & 0 & 0 \\
    1 & 1 & 0 \\
    1 & 0 & 1
\end{pmatrix}
\begin{pmatrix}
    1 & 1 & 1\\
    0 & 1 & 3\\
    0 & 2 & 8
\end{pmatrix}
\end{aligned}$$

Factorize rightmost matrix to clear the second column (below diagonal)

$$\begin{aligned}
A = 
\begin{pmatrix}
    1 & 1 & 1 \\
    1 & 2 & 4 \\
    1 & 3 & 9
\end{pmatrix}
= 
\begin{pmatrix}
    1 & 0 & 0 \\
    1 & 1 & 0 \\
    1 & 0 & 1
\end{pmatrix}
\begin{pmatrix}
    1 & 0 & 0 \\
    0 & 1 & 0 \\
    0 & 2 & 1
\end{pmatrix}
\begin{pmatrix}
    1 & 1 & 1\\
    0 & 1 & 3\\
    0 & 0 & 2
\end{pmatrix}
\end{aligned}$$

Multiply the first 2 matrices together (noting that the first column is preserved):

$$\begin{aligned}
A = 
\begin{pmatrix}
    1 & 1 & 1 \\
    1 & 2 & 4 \\
    1 & 3 & 9
\end{pmatrix}
= \qquad \;
\begin{pmatrix}
    1 & 0 & 0 \\
    1 & 1 & 0 \\
    1 & 2 & 1
\end{pmatrix} \qquad \;
\begin{pmatrix}
    1 & 1 & 1\\
    0 & 1 & 3\\
    0 & 0 & 2
\end{pmatrix}
= L U
\end{aligned}$$

And $A$ is now factored into lower triangular $L$ and upper triangular $U$.

---

More realistic example with fewer 1's:
Follow convention - Leave the pivots on diagonal in $U$ with 1's on diagonal of $L$:

$$\begin{aligned}
A = 
\begin{pmatrix}
    4 & 3 & 2 \\
    16 & 14 & 9 \\
    12 & 13 & 13
\end{pmatrix}
= 
\begin{pmatrix}
    1 & 0 & 0 \\
    4 & 1 & 0 \\
    3 & 0 & 1
\end{pmatrix}
\begin{pmatrix}
    4 & 3 & 2\\
    0 & 2 & 1\\
    0 & 4 & 7
\end{pmatrix}
\end{aligned}$$

$$\begin{aligned}
A = 
\begin{pmatrix}
    4 & 3 & 2 \\
    16 & 14 & 9 \\
    12 & 13 & 13
\end{pmatrix}
= 
\begin{pmatrix}
    1 & 0 & 0 \\
    4 & 1 & 0 \\
    3 & 0 & 1
\end{pmatrix}
\begin{pmatrix}
    1 & 0 & 0 \\
    0 & 1 & 0 \\
    0 & 2 & 1
\end{pmatrix}
\begin{pmatrix}
    4 & 3 & 2\\
    0 & 2 & 1\\
    0 & 0 & 5
\end{pmatrix}
\end{aligned}$$


$$\begin{aligned}
A = 
\begin{pmatrix}
    4 & 3 & 2 \\
    16 & 14 & 9 \\
    12 & 13 & 13
\end{pmatrix}
= \qquad \;
\begin{pmatrix}
    1 & 0 & 0 \\
    4 & 1 & 0 \\
    3 & 2 & 1
\end{pmatrix} \qquad \;
\begin{pmatrix}
    4 & 3 & 2\\
    0 & 2 & 1\\
    0 & 0 & 5
\end{pmatrix}
= L U
\end{aligned}$$

---

$L$ is lower triangular with 1's on the diagonal
$U$ is upper triangular with pivot's on main diagonal
Can be stored together in a single $n \times n$ array (in-place factorization):

- Upper triangle stores $U$ (with pivots on the diagonal)
- Off-diagonal elements of $L$ stored below the diagonal (diagonal of 1's need not be stored explicitly)
- product of pivots $= \sum_{i=0}^{n-1}a_{i,i} = Det(A)$

---

```python
# check that the factorization and det
you = np.array([[4,3,2],[0,2,1],[0,0,5]])
ell = np.array([[1,0,0],[4,1,0],[3,2,1]])
product_of_factors = np.dot(ell, you)
pivots = np.diag(you)
print("Product of factors = \n", product_of_factors, "\nU = \n", you)
print("Pivots = ", pivots)
print("Product of pivots = ", np.prod(pivots))
print("Det(A) = ", np.linalg.det(product_of_factors))
```

INSERT CORRECT OUTPUT

---

Use factorization to solve linear system as 2 triangular solves:

$A x = b \iff L U x = b \iff L y = b$ where $U x = y$

- Factor: $A \rightarrow L U$
- Solve $L y = b$
- Solve $U x = y$

Code up LU factorization and triangular solvers

---

```python
def LU_factor(A):
    m,n = A.shape
    if m != n:
        print("WARNING: Non-square input matrix")
        return
    mult = 0
    U = np.copy(A) #make a copy of the array
    #Note that U=A just makes another name for A, not a new copy of the array
    L = np.eye(n) #numpy's name for the identity matrix is "eye"
    for i in range(n): # for each row i
        for j in range(i+1,n): # for each row below row i
            mult = U[j,i]/U[i,i]
            L[j,i] = mult
            for k in range(i,n): # for each entry beyond the i^th diagonal entry
                U[j,k] = U[j,k] - mult*U[i,k] # for entries to the right, subtract multiple of term in row i         
    return L,U
```

---

Test it out:
```python
make_a(3)
L,U = LU_factor(make_a(3))
print("L= \n", L,"\nU=\n",U)
```
L= 
 [[1. 0. 0.]
 [1. 1. 0.]
 [1. 2. 1.]] 
U=
 [[1. 1. 1.]
 [0. 1. 3.]
 [0. 0. 2.]]

 ---

 ```python
# check that the factorization reproduces the input matrix
np.dot(L,U)
 ```
 array([[1., 1., 1.],
       [1., 2., 4.],
       [1., 3., 9.]])

---

```python
def upper_tri_solve(U,b):
    """
    insert docstring here
    """
    m,n = U.shape #matrix has m rows and n columns
    x=np.zeros(m) #create an array to store the solution (init to zeros)
    for i in range(m):
        row = m-i-1
        accum=0 #variable to store sum of coeffs times known entries in solution
        for j in range(i):
            accum+=U[row,j]*x[j]
        x[row]=(b[row]-accum)/U[i,i] #solve for i^th entry in solution
    return x
```

---

```python
def lower_tri_solve(L,b):
    """
    insert docstring here
    """
    m,n = L.shape #matrix has m rows and n columns
    # should really check for compatible size
    y=np.zeros(m) #create an array to store the solution (init to zeros)
    for i in range(m):
        row = i
        accum=0
        for j in range(i):           
            accum+=L[row,j]*y[j]
        y[row]=(b[row]-accum)/L[i,i] #solve for i^th entry in solution
    return y
```

---

```python
# check lower_tri_solve
b = make_b(3)
y = lower_tri_solve(L,b)
print("y=",y)
print("residual = ", np.dot(L,y)-b)
```
y= [ 1. -2.  4.]
residual =  [0. 0. 0.]

---

```python
def LU_solve(L,U,b):
    y = lower_tri_solve(L,b)
    x = upper_tri_solve(U,y)
    return x, y
```

---

```python
N = 20
A,b = make_a(N), make_b(N)
L,U = LU_factor(A)
x,y = LU_solve(L,U,b)
residual = np.dot(A,x)-b
np.set_printoptions(precision=2)
print("x=",x, "\ny=", y, "\nresidual = ", residual)
print("Norm of residual = ", np.linalg.norm(residual))
```
x= [ 1.05e+06 -3.66e+06  5.56e+06 -4.96e+06  2.94e+06 -1.24e+06  3.90e+05
 -9.36e+04  1.75e+04 -2.57e+03  2.99e+02 -2.77e+01  2.04e+00 -1.18e-01
  5.31e-03 -1.82e-04  4.58e-06 -7.98e-08  8.60e-10 -4.31e-12] 
y= [ 1.00e+00 -2.00e+00  4.00e+00 -8.00e+00  1.60e+01 -3.20e+01  6.40e+01
 -1.28e+02  2.56e+02 -5.12e+02  1.02e+03 -2.05e+03  4.10e+03 -8.19e+03
  1.64e+04 -3.28e+04  6.55e+04 -1.31e+05  2.62e+05 -5.24e+05] 
residual =  [ 1.04e-10 -2.01e-10  9.83e-09  4.74e-08 -5.19e-07 -3.04e-06 -1.61e-05
 -3.97e-05 -8.90e-05 -2.43e-04 -3.00e-04 -1.35e-03 -2.22e-03 -8.02e-03
 -3.26e-02 -3.74e-02  3.03e-02 -2.28e-02  1.47e-01  1.06e+00]

---

- Residual has entries $\approx 1$
- Norm of residual = 1.068 which is NOT small. 
- What went wrong?
- Problem becomes ill-conditioned as the matrix becomes larger. 
- A (not reliable) indicator is that $Det(A)$ becomes large.

```python
d = 1
for i in range(N):
    d *= U[i,i]
print("Det(A[", N, "]) = ", d)
```
Det(A[ 20 ]) =  5.2385873933228584e+137
Counter example ($Det(A)>>1$ but matrix solves nicely)
Return to consider reliable indicator: __condition number__

---

__Cramer's rule and Laplace expansion__
Can immediately write down the solution for any variable in the linear system $A x = b$: 
$$ x_i = \frac{Det(A_i)}{Det(A)}$$ 
- $A_i$ is $a$ with the $i^{th}$ column replaced by $b$ - - Can evaluate the determinants recursively by Laplace expansion: 
$$Det(A) = \sum_{i=0}^{n-1} (-1)^{i+j} Det(A_{i,j})$$ 
- $A_{i,j}$ is $A$ with row $i$ and column $j$ removed.
- Why not just compute the solution this way? What is the computing cost?

---

- Cramer's rule is handy for very small problems of theoretical results, BUT...
- __Using Cramer's rule (or matrix inverse) to COMPUTE the solution of a linear system is only efficient for...__

---

- Cramer's rule is handy for very small problems of theoretical results, BUT...
- __Using Cramer's rule (or matrix inverse) to COMPUTE the solution of a linear system is only efficient for..convincing someone that you do not know much about numerical methods.__
- It is expensive and ill-conditioned. 
- Whenever possible use factorization or iteration methods (coming up next)

---

## Section 2.2 - Iterative Solution Methods for $Ax = b$
There is an entire class of iterative linear solvers. 
Basic idea is a 3-step iterative scheme:

1) Make an initial estimate of the solution $x^{(0)}$.
2) Given a guess $x^{(k)}$, compute an improved estimate $x^{(k+1)}$.
3) Repeat until a stopping criterion is reached:
- Change between successive estimates is sufficiently small
- Residual error is sufficiently small
- Number of iterations reaches a specified limit

> Superscript in parentheses indicates iteration number, not an exponent.

---

__Jacobi Iteration (simplest iterative method)__
- Solve each equation/row for the corresponding solution entry 
(i.e. solve the $i^{th}$ equation for $x_i$ for `i in range(n)` )
- Update by plugging previous estimate into the "row-wise" solution.

The symbolic version involves the matrix $D$ which is the diagonal part of $A$.
In $Ax = b$, substitute $A \rightarrow D-(D-A)$ 

$\begin{aligned}
(D - (D-A))x &= b\\
Dx - (D-A)x &= b \\ 
x &= D^{-1} ((D-A)x + b)
\end{aligned}$

Evaluate RHS using previous iterate to get the Jacobi iteration update formula:
$$x^{(k+1)} = D^{-1} ((D-A)x^{(k)} + b) = D^{-1} D x^{(k)} +  D^{-1} (b - Ax^{(k)})$$

or in terms of the __residual__ $r^{(k)} = (b - A x^{(k)})$:
$$\qquad x^{(k+1)} = x^{(k)} + D^{-1} r^{(k)}$$

---

Jacobi iteration formula:  $x^{(k+1)} = x^{(k)} +  D^{-1} (b - Ax^{(k)})$

- Requires computing a matrix inverse, but __methods based on computing inverses are inefficient and ill-behaved__.
- Does that present a serious limitation for using Jacobi iteration?
- Here the matrix to invert is diagonal, so we just have to do $n$ divisions and we are OK provided that no diagonal element is close to zero.

Now let's see Jacobi iteration in action...

---

```python
# Apply Jacobi iteration to solve text Eq.(2,2,4)
A1 = np.array([[4,-1,1],[4,-8,1],[-2,1,5]])
b1 = np.array([7,-21,15])
print(A1,b1)
```
A1 =  [[ 4 -1  1]
 [ 4 -8  1]
 [-2  1  5]]
b1 =  [  7 -21  15]

```python
#Verify that the exact solution is [2,4,3]
A1DotSol = np.dot(A1, np.array([2,4,3]))
np.allclose(A1DotSol, b1), A1DotSol, b1
```
(True, array([  7, -21,  15]), array([  7, -21,  15]))

---
Write an update function based on solving line i for entry i
```python
#based on Eqs.(2.2.1)
def update_1(x,y,z):
    x_new = (7+y-z)/4.
    y_new = (-21 -4*x - z)/(-8.)
    z_new = (15+2*x-y)/5.
    return x_new,y_new,z_new
```
`update_1(1.,2.,2.)`
(1.75, 3.375, 3.0)

---

```python
iters = 6 # Reproduce some results in Table 2.1
data = np.zeros([iters,4])
x,y,z =  1., 2., 2.
for k in range(iters):
    data[k,0] = k
    data[k,1] = x
    data[k,2] = y
    data[k,3] = z
    x,y,z = update_1(x,y,z)
print(data)
```
[[0.   1.   2.   2.  ]
 [1.   1.75 3.38 3.  ]
 [2.   1.84 3.88 3.02]
 [3.   1.96 3.92 2.96]
 [4.   1.99 3.98 3.  ]
 [5.   1.99 4.   3.  ]
 [6.   2.   4.   3.  ]] (2 decimal accuracy at iteration 6; notebook $\rightarrow$ 15 places at iteration 20)

 ---

 Do things always work out so nicely?
 - Rework the same system swapping top and bottom equations
 - Results should be the same except for ordering
```python
#based on Eqs.(2.2.4) with first and last equations swapped
def update_2(x,y,z):
    x_new = (y+5*z-15)/2.
    y_new = (21+4*x+z)/8.
    z_new = y-4*x+7
    return x_new,y_new,z_new
```

---

```python
iters = 7
data = np.zeros([iters,4])
x,y,z =  1., 2., 2.
for k in range(iters):
    data[k,0] = k
    data[k,1] = x
    data[k,2] = y
    data[k,3] = z
    x,y,z = update_2(x,y,z)
print(data)
```
[[   0.      1.      2.      2.  ]
 [   1.     -1.5     3.38    5.  ]
 [   2.      6.69    2.5    16.38]
 [   3.     34.69    8.02  -17.25]
 [   4.    -46.62   17.81 -123.73]
 [   5.   -307.93  -36.15  211.28]]


 ---

 If the values in this iteration blow up, why did the first case converge nicely?
- __Diagonal dominance__
  - Magnitude of each diagonal term is larger than the sum of the magnitudes of the other coefficients in the row
  - Guarantees convergence (as in the case of `update_1`)

To deal with problems outside this narrow category, a variety of embellished algorithms have been concocted.

---

__Gauss-Seidel iteration__ (next simplest scheme)
Same as Jacobi iteration but use updated values as they become available. (See Eq. (2.2.9) for details.)
- Gauss-Seidel improves rate of convergence
- Can G-S prevent divergence in the non-diagonally-dominant case?

```python
#based on Eqs.(2.2.4) with first and last equations swapped
def update_2gs(x,y,z):
    x_new = x + (1/5.)*(15-(-2*x+y+5*z))
    y_new = y + (-1/8.)*(-21-(4*x_new-8*y+z))
    z_new = z + (1/4.)*(7-(4*x_new-y_new+z))
    return x_new,y_new,z_new
```

---

```python
# Compare with the results from update_2 where the first and last equations are swapped.
iters = 7
data = np.zeros([iters,4])
x,y,z =  1., 2., 2.
for k in range(iters):
    data[k,0] = k
    data[k,1] = x
    data[k,2] = y
    data[k,3] = z
    x,y,z = update_2gs(x,y,z)
print(data)
```
[[  0.     1.     2.     2.  ]
 [  1.     2.     3.88   2.22]
 [  2.     2.81   4.31   1.68]
 [  3.     4.38   5.03  -0.11]
 [  4.     8.24   6.73  -4.9 ]
 [  5.    18.09  11.06 -17.25]
 [  6.    43.37  22.15 -49.02]]

 ---

Does G-S iteration converge to solution?
- Not in this case; but it does diverge more slowly
- Check for diagonal dominance or see if it can be attained by permuting rows

> Continue from here with gen. derivs. & steepest descent




