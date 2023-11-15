# ABM_LSM_Optim
This repository contains codes for reconstructing dynamical systems using the least square method. 

## Overview

Suppose, we find a description of a continuous dynamical system in a form of an autonomous odrinary differential equation
$$\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x}),$$
and a number of sample points of the trajectiory $\mathbf{x}(t_i)$ is known, but the mathematical description of a derivative function $\mathbf{f}(\mathbf{x})$ is not. For example, we have a recorded  three-dimensional trajectory $\mathbf{x = (x,y,z)^\top$ as in the left pane below, shown blue. 

![Fig1](https://github.com/aikarimov/ABM_LSM_Optim/blob/main/scheme.drawio.png)

We randomly select some sample points, shown gree-yellow in the middle pane, and reconstruct equations of the system using them, obtaining:

$$\begin{cases}
\begin{aligned}
& \dot{\mathbf{x}} = -10x + 10y \\
& \dot{\mathbf{y}} = 28x - y - xz\\
& \dot{\mathbf{z}} = -2.6667z - xy\\
\end{aligned}
\end{cases}$$



## Installation
Download a zip file or via git, and then add the ABM_LSM_Optim directory to your search path:

```
>> addpath('C:\Users\...\Downloads\ABM_LSM_Optim')  
>> savepath
```
## Algorithm

The code ABM_LSM_Optim uses two basic algorithms: the least square method (LSM) for evaluating unknown coefficient of equations and the approximate Buchberger-Moller (ABM) algorithm for excluding vanishing monomials 

## Literature
The ABM and delMinorTerms routines are written following pseudocodes provided in the work

Kera, H.; Hasegawa, Y. Noise-tolerant algebraic method for reconstruction of nonlinear dynamical systems. Nonlinear Dynamics 2016, 85(1), 675-692,  https://doi.org/10.1007/s11071-016-2715-3

If you use this code or its parts in scientific work, please, cite the following papers:

1. Karimov, A.; Nepomuceno, E.G.; Tutueva, A.; Butusov, D. Algebraic Method for the Reconstruction of Partially Observed Nonlinear Systems Using Differential and Integral Embedding. Mathematics 2020, 8, 300. https://doi.org/10.3390/math8020300

2. Karimov, A.; Rybin, V.; Kopets, E.; Karimov, T.; Nepomuceno, E.; Butusov, D. Identifying empirical equations of chaotic circuit from data. Nonlinear Dyn. 2023, 111:871–886 https://doi.org/10.1007/s11071-022-07854-0

## License
MIT License
