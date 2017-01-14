#Sparse Storage Formats

This repository consists of 2 separate projects: the main project    (***Sparse Project***) and the secondary project (***NNZ vs    Dense***). 

The first project evaluates the efficiency of the new storage format for sparse matrix, called ***jnz-format***. ***jnz-format*** and the best formats of ***boost*** library ***Mapped Matrix***  and ***Compressed Matrix*** were tested in solving linear systems Ax=b with conjugate gradient (CG) method, because CG is effective for sparse matrices and integrates different storage formats with ease.  

To form test problems, 85 positively defined symmetric matrices were chosen from the sparse matrix collection of University of Florida and the same quantity of right-hand sides were generated randomly. Results confirm the superiority of ***jnz-format***. 

In the second project, the tests of 90 matrices and right-hand sides, generated for this purpose, showed that in case of 35% or more zeros in matrix,  using jnz-format instead of straightforward rectangular dense form accelerates the process of solution with CG-method.  
