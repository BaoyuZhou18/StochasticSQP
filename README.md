StochasticSQP
=============

This software is under development.  Please check back later.

Overview
--------

StochasticSQP is a software package for solving optimization problems involving a stochastic objective and deterministic constraint functions.  (Problems with deterministic objective functions can also be solved.)  It is designed to locate a minimizer of

```
min     f(x) subject to cE(x) = 0 and cI(x) <= 0, where f(x) = E[F(x,ω)]
x ∈ Rⁿ
```

where ```f```, ```cE```, and ```cI``` are continuously differentiable.  Note that the expectation could be in the form of an average of a finite number of terms.

StochasticSQP is written in Matlab.  The main author is [Frank E. Curtis](http://coral.ise.lehigh.edu/frankecurtis/).  For a list of all contributors, please see the [AUTHORS file](StochasticSQP/AUTHORS).

Citing StochasticSQP
--------------------

StochasticSQP is provided free of charge so that it might be useful to others.  Please send e-mail to [Frank E. Curtis](http://coral.ise.lehigh.edu/frankecurtis/) with success stories or other feedback.  If you use StochasticSQP in your research, then please cite the following paper:

- A. S. Berahas, F. E. Curtis, D. P. Robinson, and B. Zhou. "Sequential Quadratic Optimization for Nonlinear Equality Constrained Stochastic Optimization." [arXiv:2007.10525](https://arxiv.org/abs/2007.10525), 2020.
