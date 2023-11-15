%Reconstruction of Memristor system using ABM_LSM algorithm
%Use 3 variables for reconstruction
close all
warning off
%rng default
rng shuffle
M = 4; %dim
%simulate system
Tmax = 0.01;
h = 1e-8;
x0 = [0.0 0.10 0.10 0.0];
t = [0:h:Tmax];

opts = odeset('RelTol',1e-13,'AbsTol',1e-15);
[t,y] = ode113(@Mem3,t,x0,opts); %solve ODE
w = transpose(Mem3(0,transpose(y))); %find derivatives


Ns = length(t);


%x'=y; y'=z; z' = f(...) -> I(z) = y; I(y) = x; z' = dz
u = zeros(Ns,M);

u(:,2:4) = y(:,2:4); %assign all variables except first (it is not observable)
w(:,1) = y(:,2); %z' = x

%Reconstruct first variable through integration of second variable

val = u(:,2); %1st variable
%**** BOOL 5 ORD ********
%[iv, ind] = integrate_bool(val,0,h); %find z
%**** SIMPSON 4 ORD *****
[iv, ind] = integrate_simp(val,0,h); %find z
% %**** TRAPZ 2 ORD ****
%  ind = 1:Ns;
%  [iv] = integrate_trapz(val,0,h);
% %**** EULER 1 ORD ****
% ind = 1:Ns;
% [iv] = integrate_Euler(val,0,h);


u(ind,1) = iv; % y
yr = y(ind,:);%real

w = w(ind,:);
y = u(ind,:);
t = t(ind);
Ns = length(t);
%*************************************
%plot projection

figure(2); hold on
plot(t,y(:,1),t,y(:,2),t,y(:,3),t,yr(:,1),t,yr(:,2),t,yr(:,3));
legend('$y_1$','$y_2$','$y_3$','$\hat{y}_1$','$\hat{y}_2$','$\hat{y}_3$','interpreter','latex');
xlabel('$t$','interpreter','latex');
ylabel('$y$','interpreter','latex');

% *********** error plot ******************************
figure(3); hold on
plot(t(3:end-2),yr(3:end-2,1) - y(3:end-2,1),t(3:end-2),yr(3:end-2,2) - y(3:end-2,2),t(3:end-2),yr(3:end-2,3) - y(3:end-2,3),'-');
%plot(t,yr(:,1) - y(:,1),t(:),yr(:,2) - y(:,2),t(:),yr(:,3) - y(:,3),'-');
legend('$err_1$','$err_2$','$err_3$','interpreter','latex');
xlabel('$t$','interpreter','latex');
ylabel('$err$','interpreter','latex');

figure(1); 
plot3(y(:,1),y(:,2),y(:,3));
xlabel('\itx');
ylabel('\ity');
zlabel('\itz');

%get uniformly distributed points from the simulated attractor
sigma = deglexord(0,3,M);
[L,~] = size(sigma);

N = ceil(1.3*L); %N of data points
%N = 50;

W = zeros(N,M);
Y = zeros(N,M);

for i = 1:N %take random points from attractor
    id = ceil(rand*Ns);  %number of data point
    W(i,:) = w(id,:); %X
    Y(i,:) = y(id,:); %Y
end

%reconstruct order ideal
eps = 1e-2;
[~, O] = ApproxBM(Y, eps, sigma); %use approximate Buchberger-Moller algorithm

%Use LSM for fitting the equations with the proper coefficients
eta = 1e4;
H = cell(1,M);
T = cell(1,M);
%reconstruct each equation
errsum = 0;
for i = 1:M
    V = W(:,i);
    
    [L, ~] = size(O);
    
    E = EvalPoly(eye(L),Y,O);
    h0 = (E'*E)\(E'*V);
    %          delMinorTerms(X,V,O,eta, h0, deleteminor)
    [hi,tau] = delMinorTerms(Y,V,O,eta,h0,1); %get equation and basis
    %[hi,tau] = addMajorTerms(Y,V,O,eta); %get equation and basis
    V0 = EvalPoly(hi,Y,tau); 
    errsum = errsum + norm(V - V0); %check if norm is appropriate
    
    H{1,i} = hi;
    T{1,i} = tau;
end
errsum = errsum / M
%display equations
prettyABM(H,T)
%plot sample points
figure(1); hold on
scatter3(Y(:,1),Y(:,2),Y(:,3),23,'MarkerEdgeColor','g','MarkerFaceColor','y','LineWidth',1.5);

drawnow

%simulate results
[~,y] = ode113(@(t,x)oderecon(H,T,t,x),[0:h:Tmax],y(1,:),opts); %solve ODE

figure(1);
plot3(y(:,1),y(:,2),y(:,3),'-');




 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 