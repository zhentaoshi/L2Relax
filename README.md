# L2Relax

This is R code

* [Zhentao Shi](http://www.zhentaoshi.com/) , [Liangjun Su](http://www.mysmu.edu/faculty/ljsu/), and [Tian Xie](https://cob.sufe.edu.cn/en/Home/Teachers_Details/201?typeId=1156): [“L2 relaxation”](arxiv) (2020), Arxiv

Zhentao Shi and Zhan Gao develop the R code.

Please contact Tian Xie ([xxx](xxx)) if you have any question about the simulations and empirical applications in the paper.





### Computation Environment

R, CVXR, Rmosek

For the Matlab code, [CVX](http://cvxr.com/cvx/download/) must be installed to implement convex optimization.
[Mosek](https://www.mosek.com/resources/downloads) is recommended to facilitate CVX, but not necessary.

### Generic Functions

We add a folder `generic_functions` for the estimation procedures.
The functions are ready to take input and return output.

* `SSP_PLS_est.m` is a generic function to implement PLS.
* `PLS_example.m` is a minimum example of PLS.


