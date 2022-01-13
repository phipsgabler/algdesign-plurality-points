[![CC BY 4.0][cc-by-shield]][cc-by]

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
[GPL-2](https://cran.r-project.org/web/licenses/GPL-2).  It is a literal copy I used for reference, and just removed the parts of the code which are unneccessary for my purposes.

The [rest of the code](./src) contains a derivatve of this code, and is therefore licensed under [GPL-2](https://cran.r-project.org/web/licenses/GPL-2) as well.

The [presentation](./presentation) is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

