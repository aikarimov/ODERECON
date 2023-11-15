# ABM_LSM_Optim
This repository contains codes for reconstructing dynamical systems using the least square method. 

## Overview

Suppose, the problem is to find a description of a continuous dynamical system in a form of an autonomous odrinary differential equation:
$$\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x}).$$

Let us have a number of sample points of the trajectiory $\mathbf{x}(t_i)$, but the mathematical description of the function $\mathbf{f}(\mathbf{x})$ is unknown. 

Suppose, every line of a function $\mathbf{f}(\mathbf{x})$ is a sum of monomials with coefficients, such as

$$f_j(\mathbf{x}) = 2x + 3y^2 + 17yz \dots, $$

where $\mathbf{x} = (x,y,z, \dots) ^\top$ is a phase vector. Any combination of its entries like $x,y^2,yz$ is called a monomial, and numbers by them like $2,3,17$ are coefficients. In this case, we can use ABM_LSM_Optim to find $\mathbf{f}(\mathbf{x})$.

For example, we have a recorded three-dimensional trajectory $\mathbf{x} = (x,y,z)^\top$ as in the left pane below, shown blue. 

![Fig1](https://github.com/aikarimov/ABM_LSM_Optim/blob/main/scheme.drawio.png)

We randomly select some sample points, shown green-yellow in the middle pane, and reconstruct sparse, readable equations of the system, obtaining:

$$\begin{cases}
\begin{aligned}
& \dot{x} = -10x + 10y \\
& \dot{y} = 28x - y - xz\\
& \dot{z} = -2.6667z - xy\\
\end{aligned}
\end{cases}$$

Then, we can solve this reconstructed system with a standard matlab solver like `ode45`. The obtained trajectory is shown yellow in the right pane.

## Installation
Download a zip file or via git, and then add the ABM_LSM_Optim directory to your search path:

```matlab
>> addpath('C:\Users\...\Downloads\ABM_LSM_Optim')  
>> savepath
```
## How to use

Let us find equations of the Lorenz system. First, generate a full trajectory with a stepsize $h=0.01$ from the initial point $(0.1,0,-0.1)^\top$:
```matlab
%simulate Lorenz system
Tmax = 45;
h = 0.01;
[t,y] = ode45(@Lorenz,[0:h:Tmax],[0.1,0,-0.1]); %solve ODE
w = transpose(Lorenz(0,transpose(y))); %find derivatives
```
Then, select $N$ random points:

```matlab
%get uniformly distributed points from the simulated attractor
N = 19; %data points
M = 3; %dimension
[Ns, ~] = size(y); %get the number of data point
W = zeros(N,M); %sample derivatives
Y = zeros(N,M); %sample phase coordinates

for i = 1:N %take random points from trajectory
    id = ceil(rand*Ns);  
    W(i,:) = w(id,:); 
    Y(i,:) = y(id,:);
end
```
After that, we obtain the Lorenz equations from these 19 points using ABM_LSM_Optim. First, we use a function `PolyRegression` to obtain two cell arrays $T$ and $H$, containing all necessary information about the reconstructed system (see the section [Algorithm](https://github.com/aikarimov/ABM_LSM_Optim/tree/main#algorithm) for details):

```matlab
dmax = 2; % maximum power of the monomial
[H,T] = PolyRegression(Y,W,dmax);
```
To show the reconstruction result, we use a function `prettyABM`:

```matlab
prettyABM(H,T)
```
which outputs into the console:
```
f1 = -10*x1 + 10*x2
f2 = 28*x1 - x2 - x1*x3
f3 = -2.6667*x3 + x1*x2
```
Then, we can simulate the results using a function `oderecon`:
```matlab
[~,y] = ode45(@(t,x)oderecon(H,T,t,x),[0:h:Tmax],[0.1,0,-0.1]); %solve ODE
```
## Algorithm

First, introduce some formalism. Representation of an arbitrary $M$-dimensional function $\mathbf{f}(\mathbf{x})$ is contained in two cell arrays $H$ and $T$. Entries of $H$ are $L_i \times 1$ matrices (column vectors) of coefficients by monomials in $i$-th entry of $\mathbf{f}(\mathbf{x})$, where $L_i$ is a number of terms. Entries of $T$ are $L_i \times M$ matrices containing powers of variables in each monomial, ordered degree-lexicographically.

On the first stage of the algorithm, for each line, full matrices $T_i$ are created, containing all possible variants of powers up to $d_{max}$. Example of full degree-lexicographic ordering $\sigma$ is shown in the left of the figure, and example of the Lorenz system represented in such a way is given in the right of the figure. 

![Fig2](https://github.com/aikarimov/ABM_LSM_Optim/blob/main/handt.drawio.png)

Ordering $\sigma$ is generated with the function `deglexord(dmin,dmax,M)`. Sparse reconstruction of the equations needs eliminating all excessive terms in $T_i$ and setting correct values to entries of $H_i$. 

Suppose, we have a trajectory $Y$ which is represented as $N \times M$ matrix, with $N$ sample points and $M$ dimensions, and a derivative to the trajectory $W = \dot{Y}$ which is also represented as $N \times M$ matrix. 

First, the approximate Buchberger-Moller (ABM) algorithm runs, which excludes all monomials in $T_i$ that vanish on a given set $Y$:

```matlab
[G, O] = ApproxBM(Y, eps, sigma); %use approximate Buchberger-Moller algorithm
```
Roughly speaking, vanishing means that the monomial takes values near zero (below a threshold `eps`), and keeping it in $T_i$ makes the further step poorly conditioned. 
The function `ApproxBM` returns a border basis `G` which is not needed, and an order ideal `O` which is used as an initial guess for $T_i$.

After that, a linear regression is performed with a variant of the least square method (LSM):

```matlab
%Use LSM for fitting the equations with the proper coefficients
H = cell(1,M);
T = cell(1,M);
%reconstruct each equation
for i = 1:3
    V = W(:,i);
    [hi,tau] = delMinorTerms(Y,V,O,eta); %get equation and basis    
    H{1,i} = hi;
    T{1,i} = tau;
end
```

The function `delMinorTerms(Y,V,O,eta)` evaluates coefficients by each monomial by LSM, and then evaluates the contribution of this monomial to the whole function value on the set $Y$. If `1/N*norm(V - EvalPoly(hi,Y,tau)) <= eta `, which means the normalized error between the values of the reconstructed function and real values is not greater than `eta`, the minor term (the term which contribution is the lowest) is removed from the regression.

The function `EvalPoly(hi,Y,tau)` is used to estimate the required function described by a pair $\hi = H_i$ and $\tau = T_i$ in all points of $Y$. If we substitute the identity matrix instead of $hi$,the function will return a matrix $E$ containing values of all monomials in $tau$ in every point of $Y$:

```matlab
E = EvalPoly(eye(L),X,tau);
```

This matrix is used for estimating $\hi$ via QR decomposition:

```matlab
[Q,R] = qr(E);
Q1 = Q(:,1:L);
R1 = R(1:L,1:L);
hi = R1\(Q1'*V);
```

## Literature
The ABM and delMinorTerms routines are written following pseudocodes provided in the work

Kera, H.; Hasegawa, Y. Noise-tolerant algebraic method for reconstruction of nonlinear dynamical systems. Nonlinear Dynamics 2016, 85(1), 675-692,  https://doi.org/10.1007/s11071-016-2715-3

If you use this code or its parts in scientific work, please, cite the following papers:

1. Karimov, A.; Nepomuceno, E.G.; Tutueva, A.; Butusov, D. Algebraic Method for the Reconstruction of Partially Observed Nonlinear Systems Using Differential and Integral Embedding. Mathematics 2020, 8, 300. https://doi.org/10.3390/math8020300

2. Karimov, A.; Rybin, V.; Kopets, E.; Karimov, T.; Nepomuceno, E.; Butusov, D. Identifying empirical equations of chaotic circuit from data. Nonlinear Dyn. 2023, 111:871–886 https://doi.org/10.1007/s11071-022-07854-0

## License
MIT License
