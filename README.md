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
This funcion builds a polynomial from inputs r (even number or zero) and s (odd number) according to the following form:

![Polynomial Form](polynomial_formulas/pi_form.png?raw=true "Polynomial Form")

*fig 1: form of complete polynomial*

![Z Form](polynomial_formulas/z_form.png?raw=true "Z Form")

*fig 2: form of z terms*

![T Form](polynomial_formulas/t_form.png?raw=true "T Form")

*fig 3: form of t terms*

The function returns either only the polynomial as a Sage symbolic expression (when degree = "false", the default setting), or a list in the following format \[polynomial (as a symbolic expression), degree of the polynomial (as an integer)].

##### *find_monomials(polynomial, degree, r, s)*
This function takes a polynomial built from the form outlined by the previous function, and finds all monomials present in the expanded form of the polynomial that fit the terms described in [monomial_form.md](polynomial_formulas/monomial_form.md). The function returns a list of monomials with their coefficients as symbolic expression objects. If no monomials are found, the function returns an empty list. 

##### *sift(r, s)*
This function executes the main purpose of this module. Quite simply, it combines the previous two functions. Input is r and s as with `build_polynomial`, and the function returns a list of monomials in the same fashion as `find_monomials`. 

### Installation 
This module can be installed into your copy of Sage simply by copying partsift.py into Sage's version of Python. The location is as follows: `*/sagemath/local/lib/python2.7` with `*` being the path to where Sage is installed on your computer. Once the file is there, you can use these functions in your Sage programming at any time with `import partsift`. 
