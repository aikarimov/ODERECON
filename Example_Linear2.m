%Reconstruction of 2-dim linear system using ABM algorithm and delMinorTerms
rng default;
close all;

%simulate Linear system
a1 = -1;
a2 = -0.3;
L = [0 1; a1 a2];

x0 = [0, 1];

Tmax = 45;
h = 0.01;
[t,y] = ode45(@(t,x)(L*x),[0:h:Tmax],x0); %solve ODE
w = transpose(feval(@(x)(L*x),transpose(y))); %find derivatives

%plot x-z plane projection
figure(1); 
plot(y(:,1),y(:,2));
xlabel('\itx');
ylabel('\ity');

%get uniformly distributed points from the simulated attractor
N = 19; %data points
M = 2; %dim
[Ns, ~] = size(y);
W = zeros(N,M);
Y = zeros(N,M);

for i = 1:N %take random points from attractor
    id = ceil(rand*Ns);  %number of data point
    W(i,:) = w(id,:); %X
    Y(i,:) = y(id,:); %Y
end

%plot sample points
figure(1); hold on
scatter(Y(:,1),Y(:,2),23,'MarkerEdgeColor','g','MarkerFaceColor','y','LineWidth',1.5);

%reconstruct order ideal
eps = 1e-5;

sigma = deglexord(1,M);
[~, O] = ApproxBM(Y, eps, sigma) %use approximate Buchberger-Moller algorithm

%Use LSM for fitting the equations with the proper coefficients
eta = 1e-5;
H = cell(1,M);
T = cell(1,M);
%reconstruct each equation
for i = 1:M
    V = W(:,i);
    [hi,tau] = delMinorTerms(Y,V,O,eta); %get equation and basis
    V0 = EvalPoly(hi,Y,tau); 
    norm(V - V0) %check if norm is appropriate
    
    H{1,i} = hi;
    T{1,i} = tau;
end

%simulate results
[~,y] = ode45(@(t,x)oderecon(H,T,t,x),[0:h:Tmax],x0); %solve ODE
figure(1);
plot(y(:,1),y(:,2),'-');

%display equations
prettyABM(H,T)



 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 