# multivsos
Sums of squares decomposition of multivariate nonnegative polynomials

## Description
`multivsos` is a Maple library for computation of sums of squares (SOS) decompositions of multivariate nonnegative polynomials with rational coefficients in the (un)-constrained case. 
In the unconstrained case, `multivsos` implements a hybrid numeric-symbolic algorithm computing exact rational SOS decompositions for polynomials lying in the interior of the SOS cone. It computes an approximate SOS decomposition for a perturbation of the input polynomial with an arbitrary-precision semidefinite programming (SDP) solver. An exact SOS decomposition is obtained thanks to the perturbation terms. 
In the constrained case, `multivsos` allows to compute weighted SOS decompositions for polynomials positive over basic compact semialgebraic sets.

## Installation instructions
### Prerequisites
`multivsos` has been tested with `Maple2016` and requires: 
- the `convex` Maple package available at the address http://www-home.math.uwo.ca/~mfranz/convex
- the external `SDPA` software (SDP solver) available at the address https://sourceforge.net/projects/sdpa/files/sdpa/sdpa_7.3.8.tar.gz
- the external `SDPA-GMP` software (arbitrary-precision SDP solver) available at the address https://sourceforge.net/projects/sdpa/files/sdpa-gmp/sdpa-gmp.7.1.3.src.20150320.tar.gz


### Download
`multivsos` is maintained as a GitHub repository at the address https://github.com/magronv/multivsos.git

It can be obtained with the command:

$ git clone https://github.com/magronv/multivsos.git

### Execution and Benchmarks
From the multivsos/ directory, launch Maple and execute the following command:

`read "multivsos1.mm":`

#### Unconstrained problems

In the unconstrained case, the `multivsos1` procedure takes as input in the following order:

- the initial polynomial f 
- the perturbation magnitude (default value = 10)
- the precision of numerical Cholesky's decomposition (default value = 10)
- the number of precision digits for the SDP solver (default value = 200)
- the number of precision digits for the termination criterion of the SDP solver (default value = 30)
- the rounding precision for the output of the SDP solver (default value = 30)
- the rounding precision for the input of the SDP solver (default value = 10)
- a boolean variable for the use of arbitrary-precision SDPA-GMP solver (default value = false)
- an integer which is equal to 1 for multivsos and 2 for a previous algorithm developed by Parrilo and Peyrl in https://www.sciencedirect.com/science/article/pii/S0304397508006452 (default value = 1)

Let us consider the polynomial f := x^4 + x^3 y - 7/4 x^2 y^2 - 1/2 x y^3 + 5/2 y^4. 

To compute a sums of squares decomposition of f, you can type:

`f := x^4 + x^3*y - 7/4*x^2*y^2 - 1/2*x*y^3 + 5/2*y^4: sos,r:=multivsos1(f,2,30,200,30,2,10,false,1);`

                                                                                                            2
                         2      2             395    2      2                  2                    2      y
 sos, r := [1/12, x y - y , 0, x , 5/36, x y, ----, y , 1, x  + 1/2 x y - 4/3 y , 1, 2/3 x y + 3/4 y , 1, ----], [7]
                                              7056                                                         7


The first output sos is a list [c<sub>1</sub>,p<sub>1</sub>,...,c<sub>r</sub>,p<sub>r</sub>], where each c<sub>i</sub> is a rational number and each p<sub>i</sub> is a rational polynomial such that f admits the following weigthed SOS decomposition:

f  = c<sub>1</sub> p<sub>1</sub> <sup>2</sup> + ... + c<sub>r</sub> p<sub>r</sub> <sup>2</sup>

The second output r provides a list (with a single element in the unconstrained case) indicating the number of squares in the decomposition.

You can verify afterwards that this yields a valid nonnegativty certificate of f with the following command:

`s := 0: for i from 1 to nops(sos)/2 do s := s + sos[2*i-1]*sos[2*i]^2 od: expand (f -s);`


#### Constrained problems

In the constrained case, the `multivsos1` procedure takes as input the same parameters as in the unconstrained case as well as:

- the list of polynomial [g<sub>1</sub>,...,g<sub>m</sub>] encoding the set of inequality constraints g<sub>1</sub>(x) >= 0,...,g<sub>m</sub>(x) >= 0. For instance, the hypercube [-1,1]^2 is encoded by [1 - x^2, 1 - y^2] (default value = [])
- an integer for the maximal degree of the SOS decomposition (default value = 0)

Let us consider the polynomial f := -x^2  - 2 x y - 2 y^2 + 6 on the hypercube encoded by the list of two polynomial constraints g :=  [1 - x^2, 1 - y^2]. To compute an SOS decomposition with maximal degree of 2, you can execute the following commands:

`f:=  -x^2  - 2 * x * y - 2 * y^2  + 6: g :=  [1 - x^2, 1 - y^2]: halfdegree := 1: sos, rlist := multivsos1(f,1,10,200,30,3,100,false,1,g,halfdegree):`


                   23853407      23     130657269                              y
    sos, rlist := [---------, 1, --, x, ---------, y, 1, 1/2442, 1, x - y, 1, ----, 1, 11/7, 1, 13/7], [6, 1, 1]
                   292204836     49     291009481                             2437

The first output is a list [c<sub>0 1 </sub>, p<sub>0 1 </sub>,...,c<sub>0 r <sub>0</sub> </sub>, p<sub>0 r <sub>0</sub> </sub>,...,c<sub>m 1 </sub>, p<sub>m 1 </sub>, ...c<sub>m r <sub>m</sub> </sub>, p<sub>m r <sub>m</sub> </sub>] where each  c<sub>i j </sub> is a rational number and each p<sub>i j </sub> is a rational polynomial such that f admits the following weigthed SOS decomposition:

f  = c<sub>0 1 </sub> p  <sub>0 1 </sub> <sup>2 </sup> +...+c<sub>0 r <sub>0</sub> </sub>, p<sub>0 r <sub>0</sub> </sub> <sup>2 </sup> + g <sub> 1 </sub> [ c<sub>1 1 </sub> p<sub>1 1 </sub> <sup>2 </sup> + ... + c<sub>1 r <sub>1</sub> </sub> p<sub>1 r <sub>1</sub> </sub> <sup>2 </sup> ] +...+ g <sub> m </sub> [c<sub>m 1 </sub> p<sub>m 1 </sub> <sup>2 </sup> + ... + c<sub>m r <sub>m</sub> </sub> p<sub>m r <sub>m</sub> </sub> <sup>2 </sup>].


The second output is the list [r<sub>0 </sub>,...,r<sub>m</sub>].

You can verify afterwards that this yields a valid nonnegativty certificate of f with the following commands:

`oneg := [1,op(g)]:s := 0: idx:=0: for i from 1 to nops(g)+1 do idxi:=idx+rlist[i]:for j from idx+1 to idxi do s := s + oneg[i]*sos[2*j-1]*sos[2*j]^2: od: idx:=idxi: od: expand (f - s);`


#### Benchmarks

See the file `allbenchs.mm` 
