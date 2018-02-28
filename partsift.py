"""
Kat Cannon-MacMartin
partsift v1.5
A tool for building and polynomials and finding monomials.
for use in the paper:
'Sequences in Dihedral Groups with Distinct Partial Products'
February 20, 2018
Marlboro College
"""
from sage.calculus.var import var
from sage.combinat.partition import Partitions
from sage.calculus.functional import expand

import itertools
from sympy.utilities.iterables import multiset_permutations

def build_polynomial(r, s, degree = False):
	"""Function to build a polynomial from r and s inputs.
	Args:
	    r (int): r value.
	    s (int): s value.
	    degree (boolian): True to return degree along
			      with completed polynomial.
			      False to return only
			      polynomial.

	Returns:
	    If degree is False, returns the completed
	    polynomial variable as a sage symbolic
	    expression.
	    If degree is True, returns the list:
	    [polynomial, poly_degree] where polynomial
	    is the same as before, and poly_degree is
	    an integer.
	"""

	polynomial = 1
	"""Base for the end polynomial. Each time a new term
	is created, it wil be multiplied into this variable.
	"""

	if r % 2 == 0:
		p = r/2
		sig = 0
	else:
		p = (r-1)/2
		sig = 1

	if s % 2 == 0:
		q = (s-2)/2
		thet = 1
	else:
		q = (s-1)/2
		thet = 0
	"""p and q values are derived from r and s and are
	important mathematically for building the polynomial.
	sig and thet variables are used to indicate if r and
	s are odd or even.
	"""


	poly_degree = 0
	"""This variable is used to record the degree of the
	polynomial being built. It is most important so that
	build_polynomial() can work together nicely with
	find_monomials().
	"""

	#This block adds all basic x variables to the function.
	for j in xrange(2, r+1):
		for i in xrange(1, j):
			term = (var("x" + str(j))-var("x" + str(i)))
			polynomial = polynomial *term
			poly_degree += 1

	#This block adds all basic y variables to the function.
	for j in xrange(2, s+1):
		for i in xrange(1, j):
			term = (var("y" + str(j)) - var("y" + str(i)))
			polynomial = polynomial*term
			poly_degree += 1



	def z_adder(input):
		"""Function for use in building z and t variables.
			Args:
			    input(int): The number of the z variable
				        to be shown in terms of x.
			Returns:
			    z_added(symb exp): Series of x variables
					       added together.
		"""
		z_added = 0
		for i in xrange(1, input+1):
			z_added = z_added + var("x"+str(i))
		return z_added

	def t_adder(input, p):
		"""Function similar to z_adder(), but for t variables.
			Args:
			    input(int): The number of the t variable
					to be shown in terms of x and y.
			Returns:
			    t_added(symb exp): Series of x and y variables
					       added together.
		"""
		t_added = 0
		alty = itertools.cycle([0,1]).next
		for i in xrange(1, 2*input):
			altyvar = alty()
			if altyvar == 0:
				t_added = t_added + var("y" + str(i))
			else:
				t_added = t_added - (var("y" +str(i)))
		return t_added + z_adder(p)

	def findt(input, p, q):
		"""Function to find individual t variables in terms of x
		   and y.
			Args:
			    input(int): The number of the t variable.
			Returns:
			    t term(symb exp): The t variable expressed
					       in terms of x and y.
		"""
		if input > q+1:
			final = 0
			for l in xrange(p+1, input-q+p):
				final = final - (var("x" + str(l)))
			return t_adder(q + 1, p) + final
		else:
			return t_adder(input, p)


	def findz(input, p):
		"""Function to find individual z variables in terms of x
		   and y.
			Args:
			    input (int): The number of the z variable
				   to be shown in terms of x and y.
			Returns:
			    z term(symb exp): The z variable expressed
					      in terms of x and y.
		"""
		if input == 0:
			return 0
		else:
			if p + q + 1 > input > p:
				final = 0
				alt = itertools.cycle([0,1]).next
				for l in xrange(1, 2*(input-p) +1):
					altvar = alt()
					if altvar == 0:
						final = final + var("y" + str(l))
					else:
						final = final - (var("y" + str(l)))

				return z_adder(p) + final
			elif input == p + q + 1:
				return findt(p+q+1+sig, p, q) - var("y" + str(2*q +2))

			else:
				return z_adder(input)

	#This block adds terms with z variables to the polynomial in terms of x and y.
	for j in xrange(1, p+q+1+thet):
		for i in xrange(0, j):
			if j!=i+1:
				term = (findz(j, p)) - (findz(i, p))
				polynomial = polynomial * term
				poly_degree += 1

	#This block adds terms with t variables to the polynomial in terms of x and y.
	for j in xrange(2, p+q+2+sig):
		for i in xrange(1, j):
			if j != i+1:
				term = (findt(j, p, q)) - (findt(i, p, q))
				polynomial = polynomial * term
				poly_degree +=1

	#This block adds a missing z term in the case of s being even.
	if thet == 1:
		term = (findz(p+q+1, p) - findz(p+q, p))
		polynomial = polynomial*term
		poly_degree +=1

	if degree == True:
		result_list = [polynomial, poly_degree]
		return result_list
	else:
		return polynomial




def find_monomials(polynomial, poly_degree, r, s, form = "list", sort_type = "zero"):
	"""This function finds all monomials that, in relation to the given monomial, fit
	   the criteria of the problem.

	   ***ATTENTION: DESPITE THE VARIABLE DECLARATIONS AT THE BEGINNING OF THE FUNCTION,
	   THIS FUNCTION WILL NOT WORK CORRECTLY IF ALL VARIABLES IN THE INPUT POLYNOMIAL
	   HAVE NOT ALREADY BEEN DECLARED.***

		Args:
		    polynomial(symb exp): The polynomial to which monomials will be matched.
		    poly_degree(int): The degree of the input polynomial.
		    r(int): The r value for the input polynomial.
		    s(int): The s value for the input polynomial.
		    form(str): If form is "list", the function will return a list of monomials.
			       If form is "tup", the function will return a list of tuples with
			       the form: (monomial, coefficient).
		    sort_type(str): If sort_type is "zero", monomials are sorted by distance of
				    the coefficient from zero. If sort type is "value", monomials
				    are sorted by the coefficient's value.
		Returns:
		    monomials(list): Returns either a list of monomials or a list of tuples depending
				     on the value of form. Sorted by value or distance from zero of
				     coefficients depending on the value of sort_type. If there are
				     no working monomials, the function returns an empty list.
	"""
	#Double check that variables are declared to prevent errors when forming monomials.
	for num in xrange(1, r+1):
		var("x"+str(num))
	for num in xrange(1, s+1):
		var("y"+str(num))

	#Monomials must be checked against the expanded monomial
	polynomial_expanded = polynomial.expand()


	"""This variable will be used to create a monomial "template" with all the essential
	   variables and variables in place of exponent values.
	"""
	test_monomial = 1
	for number_of_xs in xrange(1, r+1):
		test_monomial = test_monomial*var("x"+str(number_of_xs))**var("xpower"+str(number_of_xs))
	for number_of_ys in xrange(1, s+1):
		test_monomial = test_monomial*var("y"+str(number_of_ys))**var("ypower"+str(number_of_ys))


	"""This block creates a list of placeholders for x and y variable exponents in the order
	   they appear.
	"""
	list_of_xpowers = []
	for number_of_xs in xrange(1, r+1):
		list_of_xpowers.append(var("xpower"+str(number_of_xs)))

	list_of_ypowers = []
	for number_of_ys in xrange(1, s+1):
		list_of_ypowers.append(var("ypower"+str(number_of_ys)))

	list_of_both_powers = list_of_xpowers + list_of_ypowers
	del(list_of_xpowers)
	del(list_of_ypowers)


	"""This block finds all possible exponent combinations for x and y variables
	   in a list form, then finds all possible permutations of that list.
	"""
	list_of_power_values = []

	for partition_set in Partitions(poly_degree, max_length=r+s).list():
		partition_set = list(partition_set)
		if max(partition_set)<max(r,s):
			if len(partition_set)==r+s:
				list_of_power_values.append(partition_set)
			else:
				while len(partition_set)<r+s:
					partition_set.append(0)
				list_of_power_values.append(partition_set)


	"""Now, every list of exponent possibilites is subbed into the monomial template,
	   then checked against the expanded polynomial. If it is present in the expanded
	   form, it is added to the final list.
	"""
	workable_monomials = []
	for power_values in list_of_power_values:
		for power_set in multiset_permutations(power_values):
			temp_test_monomial = test_monomial
			power_set=list(power_set)
			if max(power_set[0:r])<r and max(power_set[r:r+s])<s:
				temp_test_monomial = temp_test_monomial.subs({list_of_both_powers[i]:power_set[i] for i in xrange(0, r+s)})
				if polynomial_expanded.coefficient(temp_test_monomial) not in [0]:
					workable_monomials.append((polynomial_expanded.coefficient(temp_test_monomial)*temp_test_monomial, int(polynomial_expanded.coefficient(temp_test_monomial))))

	#Changes output depending on prefernces set at the outset.
	if len(workable_monomials) > 0:
		if sort_type == "zero":
			workable_monomials = sorted(workable_monomials, key = lambda tup: abs(tup[1]))
		else:
			workable_monomials = sorted(workable_monomials, key = lambda tup: tup[1])

	if  form == "list":
		for i in xrange(0, len(workable_monomials)):
			workable_monomials[i] = workable_monomials[i][0]

	return workable_monomials


def sift(r, s, form = "list", sort_type = "zero"):
	"""This function builds a polynomial according to r and s values
	   then finds all monomials present in the expanded form of the
	   polynomial that fit the criteria of the problem.
		Args:
		    r(int): The r value.
		    s(int): The s value.
		    form(str): Assign as "list" for the function to return
			       a list of monomials. Assign as "tup" for the
			       list to contain tuples with the form:
			       (monomial, coefficient of monomial).
		    sort_type(str): Assign as "zero" to sort monomial output
				    by distance of the coefficient from zero.
				    assign as "value" to sort monomial output
				    by the value of coefficients.
		Returns:
		    List of monomials(list): List of working monomials(symb exp)
					     or tuples containing (monomial, coefficient)
					     depending on the value of form. If no monomials
					     fit the criteria, an empty list is returned.
	"""

	poly_input = build_polynomial(r, s, degree = True)
	return find_monomials(poly_input[0], poly_input[1], r, s, form = form, sort_type = sort_type)



