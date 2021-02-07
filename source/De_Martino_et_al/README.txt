Supporting material to the article 
STATISTICAL MECHANICS FOR METABOLIC NETWORKS IN STEADY-STATE GROWTH


boundscoli.dat
Flux Bounds of the E. coli catabolic core model iAF1260 in a glucose limited minimal medium. 


polcoli.dat
Matrix enconding the polytope of the E. coli catabolic core model iAF1260 in a glucose limited minimal medium, 
obtained from the soichiometric matrix by standard linear algebra  (reduced row echelon form).

ellis.dat
Approximate Lowner-John ellipsoid rounding the polytope of the E. coli catabolic core model iAF1260 in a glucose limited minimal medium
obtained with the Lovasz method.

point0.dat
Center of the approximate Lowner-John ellipsoid rounding the polytope of the E. coli catabolic core model iAF1260 in a glucose limited minimal medium
obtained with the Lovasz method.



lovasz.cpp  
This c++ code file receives in input the polytope of the feasible steady states of a metabolic network, 
(matrix and bounds), and it gives in output an approximate Lowner-John ellipsoid rounding the polytope
with the Lovasz method 
NB inputs are referred by defaults to the catabolic core of the E.Coli network iAF1260. 
For further details we refer to  PLoS ONE 10.4 e0122670 (2015).

sampleHRnew.cpp  
This c++ code file receives in input the polytope of the feasible steady states of a metabolic network, 
(matrix and bounds), the ellipsoid rounding the polytope, a point inside and  
it gives in output a max entropy sampling at fixed average growth rate 
of the steady states by performing an Hit-and-Run Monte Carlo Markov chain.
NB inputs are referred by defaults to the catabolic core of the E.Coli network iAF1260. 
For further details we refer to  PLoS ONE 10.4 e0122670 (2015).






