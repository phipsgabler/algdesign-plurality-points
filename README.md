# Plurality Points

This repository contains the [slides](./presentation/presentation.pdf) and some code for a seminar talk I gave 
about the paper [1]. 

## Bibliography

These also include the most important references for understanding the main paper.

1. Y.-W. Wu, W.-Y. Lin, H.-L. Wang, and K.-M. Chao, “Computing Plurality Points and Condorcet Points in Euclidean Space,” in International Symposium on Algorithms and Computation, 2013, pp. 688–698.
2. P. J. Rousseeuw and A. Struyf, “Computing location depth and regression depth in higher dimensions,” Statistics and Computing, vol. 8, no. 3, pp. 193–203, 1998.
3. N. Megiddo, “Linear programming in linear time when the dimension is fixed,” Journal of the ACM (JACM), vol. 31, no. 1, pp. 114–127, 1984.
4. T. H. Cormen, C. E. Leiserson, R. L. Rivest, and C. Stein, Introduction to Algorithms, 3rd ed. 2009.


## License

The file [`depth.f`](./depth.f) contains Fortran code from the R implementation of the algorithm of Rousseeuw and Struyf
(see [here](https://CRAN.R-project.org/package=depth)), which is licensed under
[GPL-2](https://cran.r-project.org/web/licenses/GPL-2).  I removed the parts of the code which are unneccessary for 
my purposes.

My code is licensed under the MIT license.  With the slides, do anything you want.

