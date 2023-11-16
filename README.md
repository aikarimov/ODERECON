# ODERECON
This repository contains codes for reconstructing dynamical systems using the least square method. 

## Overview

Suppose, the problem is to find a description of a continuous dynamical system in a form of an autonomous odrinary differential equation:
$$\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x}).$$

Let us have a number of sample points of the trajectiory $\mathbf{x}(t_i)$, but the mathematical description of the function $\mathbf{f}(\mathbf{x})$ is unknown. 

Suppose, every line of a function $\mathbf{f}(\mathbf{x})$ is a sum of monomials with coefficients, such as

$$f_j(\mathbf{x}) = 2x + 3y^2 + 17yz \dots, $$

where $\mathbf{x} = (x,y,z, \dots) ^\top$ is a phase vector. Any combination of its entries like $x,y^2,yz$ is called a monomial, and numbers by them like $2,3,17$ are coefficients. In this case, we can use ODERECON to find $\mathbf{f}(\mathbf{x})$.

For example, we have a recorded three-dimensional trajectory $\mathbf{x} = (x,y,z)^\top$ as in the left pane below, shown blue. 

![Fig1](https://github.com/aikarimov/ODERECON/blob/main/scheme.drawio.png)

We randomly select some sample points, shown green-yellow in the middle pane, and reconstruct sparse, readable equations of the system:

$$\begin{cases}
\begin{aligned}
& \dot{x} = -10x + 10y \\
& \dot{y} = 28x - y - xz\\
& \dot{z} = -2.6667z - xy\\
\end{aligned}
\end{cases}$$

Then, we can solve this reconstructed system with a standard matlab solver like `ode45`. Such a solution is shown yellow in the right pane.

## Installation
Download a zip file or via git, and then add the ODERECON directory to your search path:

```matlab
>> addpath('C:\Users\...\Downloads\ODERECON')  
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
After that, we can obtain the Lorenz equations from these $N$ toy data points using ODERECON. First, we use a function `PolyRegression` to obtain two cell arrays $T$ and $H$, containing all necessary information about the reconstructed system (see the section [Algorithm](https://github.com/aikarimov/ODERECON/tree/main#algorithm) for details):

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

![Fig2](https://github.com/aikarimov/ODERECON/blob/main/handt.drawio.png)

Ordering $\sigma$ is generated with the function `deglexord(dmin,dmax,M)`. Sparse reconstruction of the equations needs eliminating all excessive terms in $T_i$ and setting correct values to entries of $H_i$. 

Suppose, we have a trajectory $Y$ which is represented as $N \times M$ matrix, with $N$ sample points and $M$ dimensions, and a derivative to the trajectory $W = \dot{Y}$ which is also represented as $N \times M$ matrix. 

First, the approximate Buchberger-Moller (ABM) algorithm runs, which excludes all monomials in $T_i$ that vanish on a given set $Y$:

```matlab
[~, O] = ApproxBM(Y, eps, sigma); %use approximate Buchberger-Moller algorithm
```
Roughly speaking, vanishing means that the monomial takes values near zero on the entire set $Y$, and keeping it in $T_i$ makes the problem of finding $H_i$ poorly conditioned. The function `ApproxBM` returns a border basis as the first unused argument, and an order ideal `O`. The latter is utilized as an initial guess for $T_i$. 

Then, a time comes for a simple trick to find $H_i$. In many parts of a code, a function `EvalPoly(chi,Y,tau)` is used to estimate a value returned by the function described with a pair $\chi = H_i$ and $\tau = T_i$ in a point, or a set of points $Y$. If we substitute the identity matrix instead of $\chi$,the function `EvalPoly` returns a matrix $E$ containing values of all monomials in $\tau$ in every point of $Y$:

```matlab
E = EvalPoly(eye(L),X,tau);
```

This matrix is used for estimating $\chi$ via LSM using QR decomposition:

```matlab
[Q,R] = qr(E);
Q1 = Q(:,1:L);
R1 = R(1:L,1:L);
chi = R1\(Q1'*V);
```

If we do not need a sparse regression, the code stops its work. Otherwise, a function `delMinorTerms` runs for each dimension:

```matlab
%Use LSM for fitting the equations with the proper coefficients
H = cell(1,M);
T = cell(1,M);
%reconstruct each equation
for i = 1:M
    V = W(:,i);
    [chi,tau] = delMinorTerms(Y,V,O,eta); %get equation and basis    
    H{1,i} = chi;
    T{1,i} = tau;
end
```

The function `delMinorTerms(Y,V,O,eta)` estimates coefficients by each monomial as shown before, and then evaluates the contribution of this monomial to the whole function on the set $Y$. While `1/N*norm(V - EvalPoly(chi,Y,tau)) <= eta `, i.e. the normalized error between the values of the reconstructed function and real values is not greater than `eta`, the current term which contribution is the lowest is removed from the regression, new $\chi, \tau$ are found, and the procedure is repeated. In the end, the regression becomes sparce.

Optionally, iteratively reweighted least squares (IRLS) method can be used instead of an ordinary LSM, or the LASSO regression, which in some cases gives more sparse solution.

The described algorithm is implemented in the function `PolyRegression`. 

## Examples
A number of examples with use cases is provided
### Example_FHN
This example illustrates reconstructing FitzHugh-Nagumo system from 20 random data points. Running code as is, we obtain the following problem solution:
```
f1 =  x2
f2 = - x1 + 4*x2 - x2^3
```
This result is obtained when a derivative point set $W$ is obtained analytically:
```matlab
w = transpose(FHN(0,transpose(y)));
```
You can change 1 to 0 in the following code line:
```matlab
derivativecalc = 1; % Set 1 to find derivatives analytically, set 0 to find derivatives numerically
```
After that, the derivative will be found numerically with a 4-th order finite difference derivative:
```matlab
w = diff4(y,t); 
```
The result will lose its accuracy showing high sensitivity of the algorithm to numerical and truncation errors:
```
f1 = 0.99996*x2
f2 = -0.99811*x1 + 3.9906*x2 - 0.99748*x2^3
```
Adding strict accuracy options for the ode solver and replacing the `ode45` solver to the more accurate `ode78` or `ode113` improves this result, please, check it out.
### Example_Linear2

This code illustrates reconstruction of a stable 2-dimensional linear system with exponential decay. Derivatives are found analytically. The output of the code is
```
f1 =  x2
f2 = - x1 - 0.3*x2
```

### Example_Lorenz

Reconstructing the classical Lorenz attractor with a workflow shown earlier. The result of the code is
```
f1 = -10*x1 + 10*x2
f2 = 28*x1 - x2 - x1*x3
f3 = -2.6667*x3 + x1*x2
```

### Example_Lorenz_1var

Reconstructing the Lorenz attractor given only one state variable $z$. For better results, `ode113` solver is used, and additional options for tolerances are added:
```matlab
opts = odeset('RelTol',1e-13,'AbsTol',1e-15);
[t,y] = ode113(@Lorenz2,[0:h:Tmax],[0.1,0,-0.1],opts); %solve ODE
```
Once we have only $z$ variable, we must reconstruct two other variables. Usually, a direct application of Takens's theorem is performed, and the system is reconstructed in delay variables. Alternatively, derivatives are also often used. Here, we reconstruct two missing variables via numerical integration of $z$ using Bool's rule:
```matlab
val = u(:,3); %3rd variable
[iv, ind] = integrate_bool(val,0,h);
[iv2, ind2] = integrate_bool(iv,0.1,4*h);
```
We can find analytically that the resulting formula should contain negative powers, so we generate a degree-lexicographic ordering starting from $-1$:
```matlab
sigma = deglexord(-1,3,3);
```
Then, we do not use `PolyRegression` but explicitly apply `delMinorTerms` to each equation to find $H$ and $T$. This is done just for illustration, and all that bulcky code can be rewritten with `PolyRegression`. The result of the code is
```
f1 =  x2
f2 =  x3
f3 = 720*x1 - 29.3333*x2 - 13.6667*x3 + 11*x1^-1*x2^2 + x1^-1*x2*x3 - 10*x1^3 - x1^2*x2
```
For more detail, we refer to the [original publication](https://doi.org/10.3390/math8020300).

### Example_Mem3var

Reconstructing B. Muthuswamy’s circuit equations from a three-dimensional data series, with one missing state variable. The result of the code is
```
f1 =  x2
f2 = 24519*x2 + 73530*x3 - 0.086055*x1^2*x2
f3 = 7353*x2 - 7353*x3 - 14700000*x4
f4 = 55.55*x3
```
For more detail, we also refer to the [original publication](https://doi.org/10.3390/math8020300).

### Example_RingTest

This example shows how ODERECON deals with a conservative system. The output of the code is
```
f1 = -9.9865e-05*x2
f2 = 1.0006*x1
```
The inaccurate coefficients are due to numerical differentiation `dX = diff(X,1,1)`. Replacing it with the following code:
```matlab
dX = diff4(X,tspan)
```
results in accurate equations:
```
f1 = -0.0001*x2
f2 =  x1
```
## Literature
The `ApproxBM` and `delMinorTerms` functions are written following pseudocodes provided in the work

Kera, H.; Hasegawa, Y. Noise-tolerant algebraic method for reconstruction of nonlinear dynamical systems. Nonlinear Dynamics 2016, 85(1), 675-692,  https://doi.org/10.1007/s11071-016-2715-3

If you use this code or its parts in scientific work, please, cite the following papers:

1. Karimov, A.; Nepomuceno, E.G.; Tutueva, A.; Butusov, D. Algebraic Method for the Reconstruction of Partially Observed Nonlinear Systems Using Differential and Integral Embedding. Mathematics 2020, 8, 300. https://doi.org/10.3390/math8020300

2. Karimov, A.; Rybin, V.; Kopets, E.; Karimov, T.; Nepomuceno, E.; Butusov, D. Identifying empirical equations of chaotic circuit from data. Nonlinear Dyn. 2023, 111:871–886 https://doi.org/10.1007/s11071-022-07854-0

## License
MIT License
