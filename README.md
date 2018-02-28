# partsift
This is a simple set of Python 2.7 functions designed for use in the SageMath prompt or SageMath programs.
These functions were designed to aid in testing theories for the research paper:
            
            Sequences in Dihedral Groups with Distinct Partial Products
                   
                  M. A. Ollis, Marlboro College
### Required packages
partsift.py requires elements of the following python and sage packages to work:
- itertools
- sympy
- sage.calculus
- sage.combinat

### Explanation of functions:

##### *build_polynomial(r, s, degree = "false")*
This funcion builds a polynomial from inputs r (even number or zero) and s (odd number) according to the form described in [this PDF](https://drive.google.com/file/d/0B4oI_aO3gvmWRnVvMzlRa2JRYk14TGtlODhqb1hqSDVJUUx3/view?usp=sharing).

The function returns either only the polynomial as a Sage symbolic expression (when degree = "false", the default setting), or a list in the following format \[polynomial (as a symbolic expression), degree of the polynomial (as an integer)].

##### *find_monomials(polynomial, degree, r, s, form = "list", sort_type = "zero")*
This function takes a polynomial built from the form outlined by the previous function, and finds all monomials present in the expanded form of the polynomial that fit the terms described in [monomial_form.md](polynomial_formulas/monomial_form.md). The function has a variety of options for returns, chosen by the `form` and `sort_type` variables. By default, it returns a list of monomials ordered by the distance of their coefficient from zero. Setting `form = "tup"` or `form = "tuple"` will change the output to a list of tuples containing `(monomial*coefficient, coefficient)` and setting `sort_type = "value"` will order the list by the values of coefficients. 

##### *sift(r, s, form = "list", sort_type = "zero")*
This function executes the main purpose of this module. Quite simply, it combines the previous two functions. Input is r and s as with `build_polynomial`, and the function returns a list of monomials in the same fashion as `find_monomials`. `sift` also takes `form` and `sort_type` arguments in the same manner as `find_monomials`. 

### Installation 
This module can be installed into your copy of Sage simply by copying partsift.py into Sage's version of Python. The location is as follows: `*/sagemath/local/lib/python2.7` with `*` being the path to where Sage is installed on your computer. Once the file is there, you can use these functions in your Sage programming at any time with `import partsift`. 
