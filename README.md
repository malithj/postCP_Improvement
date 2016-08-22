# postCP_Improvement
postCP Package improvement for GSOC 2016

The project aimed at improving the postCP package and making it available on CRAN again. The snag that prevented the package from being updated is the recent requirement that in the R code, .C() calls require DUP=TRUE arguments, and .Call() is suggested instead of .C(). The implementation of postCP package required that it's done in the most R compliant way. Separating out the model specific implementation and the core implementation. The core part is implemented in C++ for speed in calculations. The project page can be found here. https://github.com/rstats-gsoc/gsoc2016/wiki/postCP-change-point-detection
