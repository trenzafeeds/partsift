"""
Kat Cannon-MacMartin
partsift v1.1
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

def build_polynomial(r, s, degree = "false"):
	polynomial = 1

	p = r/2
	q = (s-1)/2
	poly_degree = 0


	for j in xrange(2, r+1):
		for i in xrange(1, j):
			term = (var("x" + str(j))-var("x" + str(i)))
			polynomial = polynomial *term
			poly_degree += 1


	for j in xrange(2, s+1):
		for i in xrange(1, j):
			term = (var("y" + str(j)) - var("y" + str(i)))
			polynomial = polynomial*term
			poly_degree += 1



	def z_adder(input):
		z_added = 0
		for i in xrange(1, input+1):
			z_added = z_added + var("x"+str(i))
		return z_added


	def findz(input, p):
		if input == 0:
			return 0
		else:
			if input > p:
				final = 0
				alt = itertools.cycle([0,1]).next
				for l in xrange(1, 2*(input-p) +1):
					altvar = alt()
					if altvar == 0:
						final = final + var("y" + str(l))
					else:
						final = final - (var("y" + str(l)))

				return z_adder(p) + final
			else:
				return z_adder(input)


	for j in xrange(1, p+q+1):
		for i in xrange(0, j):
			if j!=i+1:
				term = (findz(j, p)) - (findz(i, p))
				polynomial = polynomial * term
				poly_degree += 1


	def t_adder(input, p):
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
		if input > q+1:
			final = 0
			for l in xrange(p+1, input-q+p):
				final = final - (var("x" + str(l)))
			return t_adder(q + 1, p) + final
		else:
			return t_adder(input, p)


	for j in xrange(2, p+q+2):
		for i in xrange(1, j):
			if j != i+1:
				term = (findt(j, p, q)) - (findt(i, p, q))
				polynomial = polynomial * term
				poly_degree +=1

	if degree != "false" and degree != "False":
		result_list = [polynomial, poly_degree]
		return result_list
	else:
		return polynomial




def find_monomials(polynomial, poly_degree, r, s, form = "list", sort_type = "zero"):
	for num in xrange(1, r+1):
		var("x"+str(num))

	for num in xrange(1, s+1):
		var("y"+str(num))

	polynomial_expanded = polynomial.expand()

	test_monomial = 1

	for number_of_xs in xrange(1, r+1):
		test_monomial = test_monomial*var("x"+str(number_of_xs))**var("xpower"+str(number_of_xs))
	for number_of_ys in xrange(1, s+1):
		test_monomial = test_monomial*var("y"+str(number_of_ys))**var("ypower"+str(number_of_ys))

	list_of_xpowers = []
	for number_of_xs in xrange(1, r+1):
		list_of_xpowers.append(var("xpower"+str(number_of_xs)))

	list_of_ypowers = []
	for number_of_ys in xrange(1, s+1):
		list_of_ypowers.append(var("ypower"+str(number_of_ys)))

	list_of_both_powers = list_of_xpowers + list_of_ypowers
	del(list_of_xpowers)
	del(list_of_ypowers)

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

	workable_monomials = []
	for power_values in list_of_power_values:
		for power_set in multiset_permutations(power_values):
			temp_test_monomial = test_monomial
			power_set=list(power_set)
			if max(power_set[0:r])<r and max(power_set[r:r+s])<s:
				temp_test_monomial = temp_test_monomial.subs({list_of_both_powers[i]:power_set[i] for i in xrange(0, r+s)})
				if polynomial_expanded.coefficient(temp_test_monomial) not in [0]:
					workable_monomials.append((polynomial_expanded.coefficient(temp_test_monomial)*temp_test_monomial, int(polynomial_expanded.coefficient(temp_test_monomial))))

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
	poly_input = build_polynomial(r, s, degree = "true")
	return find_monomials(poly_input[0], poly_input[1], r, s, form = form, sort_type = sort_type)



