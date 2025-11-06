
import numpy as np
import matplotlib.pyplot as plt
import itertools



## problem 1
def digits_to_num(sign, exponent, digits, beta=2):
	if sign == 0:
		s = 1
	else:
		s = -1
	sum = 0
	for i in range(len(digits)):
		sum= pow(1/2,i)*digits[i]+sum
	
	return(s*sum*pow(beta,exponent))

# print(digits_to_num(1, 3, np.array([0, 0,0,1]), beta = 2))
# print(np.array([0,0,0,1]))

## problem 2
## problem 2a
def mini_gamut():
    sign = 0
    e = 0
    bata = 2
    d = itertools.product([0, 1], repeat=4)
    digits = [i for i in d]
    save_list = []
    for x in range(0,pow(2,4)):
        save_list.append(digits_to_num(sign, e, digits[x],bata))
    return sorted(save_list)
# print(mini_gamut())

## problem 2b

def positive_gamut(e_min, e_max):
    sign = 0
    bata = 2
    d = itertools.product([0, 1], repeat=4)
    digits = [i for i in d]
    save_list = []
    for e in range(e_min,e_max+1):
        for x in range(0,pow(2,4)):
            save_list.append(digits_to_num(sign, e, digits[x],bata))
    return sorted(save_list)
# print(positive_gamut(-4,3))


## problem 2c

def gamut(e_min, e_max):
    bata = 2
    d = itertools.product([0, 1], repeat=4)
    digits = [i for i in d]
    save_list = []
    for sign in range(0,2):
        for e in range(e_min,e_max+1):
            for x in range(0,pow(2,4)):
                save_list.append(digits_to_num(sign, e, digits[x],bata))
    return sorted(save_list)
# print(gamut(-4,3))


## problem 3
def round_value(val, gamut):
    a_list= gamut
    b_list=[]
    for i in range(0,len(gamut)):
        b_list.append(abs(a_list[i]-val))
    min_index = np.argmin(b_list)            
    return a_list[min_index]
# val = 2
# gamut=gamut(-4,3)
# print(round_value(val, gamut))


## problem 4a

def gamut_linspace(count, gamut):
    x = np.linspace(gamut[0], gamut[-1], count)
    return x
# count = 5
# print(gamut_linspace(count,gamut))

##problem 4b

def round_values(vals, gamut):
    save_list=[]
    for k in range(0,len(vals)):
        save_list.append(round_value(vals[k], gamut))
    return save_list
# vals=[7,5.3,-4.99]
# print(round_values(vals, gamut))


# round_values(vals, gamut)

# vals = np.array(sorted(gamut_linspace(1000, np.array(gamut))))

def plot_rounded_error(vals, rounded_vals):
	'''
	args:
		vals: array of float values
		rounded_vals: array of float values rounded to the closest representable number

	return:
	'''

	# #uncomment when `e_abs` and `e_rel` are implemented
	e_abs = abs(rounded_vals - vals)
	e_rel = abs((rounded_vals-vals)/vals)

	fig, ax = plt.subplots(figsize = (12, 6))

	ax.plot(vals, e_abs, color = 'blue', label = 'Absolute error')
	ax.plot(vals, e_rel, color = 'black', label = 'Relative error')
	plt.legend()
	plt.ylim([0,0.6])
	plt.show()
        
# vals = np.array(sorted(gamut_linspace(1000, gamut(-4,3))))

# rounded_vals=round_values(vals, gamut(-4,3))
# plot_rounded_error(vals,rounded_vals)



# def p1():
# 	""" 
# 	1) Consider a 1-byte floating-point system with 1 sign bit, 4 precision bits (for the mantissa), and 3 exponent bits corresponding to:
# 		beta = 2, p = 4, e in {-4,-3,-2,-1,0,1,2,3}.
# 	Implement `digits_to_num`. It should return the represented number corresponding to the specific set of input values.
# 	The function calls below should produce the largest and smallest numbers the system can represent. What should these values be? Check
# 	that your function produces the appropriate values. 
# 	"""

	# sign = 0

	# e = -4
	# digits = [1,1,1,1]
	# print("largest: ", digits_to_num(sign, e, digits))

	# e = 3
	# digits = [1,0,0,0]
	# print("smallest: ", digits_to_num(sign, e, digits))

# def p2():
# 	"""
# 	2a) Implement `mini_gamut()`. It should return an array of the representable positive numbers for e = 0 and beta = 2.
# 	"""

# 	"""
# 	2b) Implement `positive_gamut(e_min, e_max)`. It should return a sorted array of all positive representable values. Use e_min and e_max
# 	from problem 1.
# 	"""

# 	"""
# 	2c) Implement `gamut(e_min, e_max)`. It shoud return a sorted array of all representable values for the example system from problem 1. Use
# 	e_min and e_max from problem 1.
# 	"""

# 	pass

# def p3():
# 	"""
# 	3) Implement `round_value(val, gamut)`. It should return the closest value from `gamut` to the numerical input value, `val`. Use the gamut generated
# 	from problem 2.
# 	"""

# 	pass

# def p4():
# 	"""
# 	4a) Implement `gamut_linspace(count, gamut)`. It should return an array with `count` values linearly spaced between the min and max gamut value.
# 	"""

# 	"""
# 	4b) Implement `round_values(vals, gamut)`. It should return an array of values which have been rounded.
# 	"""

# 	"""
# 	4c) Use `plot_rounded_error(vals, rounded_vals)` along with the functions from 4a and 4b to generate an array of 1000 values and 1000 rounded values
# 	and plot the relative and absolute errors between the rounded and unrounded results.  `plot_rounded_error` is mostly implemented, you only need to
# 	implement the error computations. Use the gamut generated from problem 2.
# 	"""

# 	pass

def p5():
	"""
	5) What is the value of machine epsilon for the floating-point system from problem 1? Explain how it plays a prominent role in one of the error plots.
	"""

	return """The machine Epsilon of N bits data type is 2^(-N). Our system is 4 bits data type, so the machine epsilon is 2^(-4).
The machine Epsilon is the maximum and an upper bound on the relative error of the system.
So, the machine epsilon is useful to determine the convergence of many iterative numerical algorithms."""
print(p5())



def p6():
	"""
	6) Do Problem 3.4 in Schorghofer
	"""

	return """(6.a)In the function f(x)= x-sqrt(x), there is sensitivity when x is close to 1.
Therefore, to keep away from loss of precision, we should take another function to replace sqrt(x). Choose the Taylor series expansion of sqrt(x) at point x=1. 

(6.b)In the function f(x)= (cos(x))^2-(sin(x))^2, there is sensitivity around the point at x=~(n+1/2)pi and n is an integer from 0. 
Therefore, to avoid this problem, we should take another function that equals (cos(x))^2-(sin(x))^2, and this alternative function is cos(2x)."""
print(p6())
	

# if __name__ == "__main__":
# 	p1()
# # 	p2()
# # 	p3()
# # 	p4()
# 	p5()
# 	p6()
