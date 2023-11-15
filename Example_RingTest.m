%Reconstruct ringtest

close all

A = sin([0.5*pi/2:0.01:10]'); %x value
iA = cos([0.5*pi/2:0.01:10]'); %integral
tspan = 0:1:length(A)-1;

h = 1;
Tmax = 1000;

figure(2);
plot(tspan,A);

tstart = 1;
p = [0.01; 1];

%construct data array
X = [p(1)*iA(tstart:end-1), p(2)*A(tstart:end-1)];

%build attractor
figure(1);
plot(X(:,1),X(:,2));
xh = xlabel('\itx');
yh = ylabel('\ity');
zh = zlabel('\itz');


dX = diff(X,1,1);
X = X(1:end-1,:); %shorten the data array to be equal to the length of dX 

%get uniformly distributed points
N = 10; %data points
M = 2; %dim
[Ns, ~] = size(X);
W = zeros(N,M);
Y = zeros(N,M);

for i = 1:N
    id = ceil(rand*Ns); %number of data point
    W(i,:) = dX(id,:); %dX, target vector
    Y(i,:) = X(id,:); %X
end

%plot points
figure(1); hold on
scatter(Y(:,1),Y(:,2),23,'MarkerEdgeColor','g','MarkerFaceColor','y','LineWidth',1.5);

%reconstruct order ideal
eps = 0.01;
sigma = [0 0; 1 0; 0 1; 2 0; 1 1; 0 2; 3 0 ; 2 1;  1 2;  0 3];
[~, O] = ApproxBM(Y, eps, sigma)

eta = 0.01;
H = cell(1,M);
T = cell(1,M);
%reconstruct each equation
for i = 1:M
    V = W(:,i);
    eps = 0.01;
    [hi,tau] = delMinorTerms(Y,V,O,eta); %get equation and basis
    V0 = EvalPoly(hi,Y,tau); 
    norm(V - V0) %check if norm is appropriate
    
    H{1,i} = hi;
    T{1,i} = tau;
end

[~,y] = ode45(@(t,x)oderecon(H,T,t,x),[0:h:Tmax],X(1,:)); %solve ODE
figure(1);
plot(y(:,1),y(:,2),'.-');

%display equations
prettyABM(H,T)