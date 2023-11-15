%Reconstruction of Lorenz system using ABM algorithm and delMinorTerms
%Use only 1 variable for reconstruction
close all;
rng shuffle
M = 3; %dim
%simulate Lorenz system
Tmax = 45;
h = 1e-4;
opts = odeset('RelTol',1e-13,'AbsTol',1e-15);
[t,y] = ode113(@Lorenz2,[0:h:Tmax],[0.1,0,-0.1],opts); %solve ODE
%[t,y] = ode45(@Lorenz2,[0:h:Tmax],[0.1,0,-0.1],opts); %solve ODE
w = transpose(Lorenz2(0,transpose(y))); %find derivatives

Ns = length(t);

%find using integration
%*************************************
%x'=y; y'=z; z' = f(...) -> I(z) = y; I(y) = x; z' = dz
u = zeros(Ns,M);
w = zeros(Ns,M);

u(:,3) = y(:,3); %z
w(:,3) = diff4(y(:,3))/h; % z'

val = u(:,3); %3rd variable

%**** MILNE 4 ORD *****
[iv, ind] = integrate_bool(val,0,h);
[iv2, ind2] = integrate_bool(iv,0.1,4*h);
%**** SIMPSON 3 ORD *****
% [iv, ind] = integrate_simp(val,0,h);
% [iv2, ind2] = integrate_simp(iv,0.1,2*h);
% %**** TRAPZ 2 ORD ****
% ind = 1:Ns; ind2 = 1:Ns;
% [iv] = integrate_trapz(val,0,h);
% [iv2] = integrate_trapz(iv,0.1,h);
% %**** EULER 1 ORD ****
% ind = 1:Ns; ind2 = 1:Ns;
% [iv] = integrate_Euler(val,0,h);
% [iv2] = integrate_Euler(iv,0.1,h); 

u(ind(ind2),2) = iv(ind2); % y
w(ind(ind2),2) = u(ind(ind2),3); %y'=z

u(ind(ind2),1) = iv2; % x
w(ind(ind2),1) = u(ind(ind2),2); %x'=y

yr = y(ind(ind2),:);%real

w = w(ind(ind2),:);
y = u(ind(ind2),:);
t = t(ind(ind2));
%*************************************
%plot projection

figure(2); hold on
plot(t,y(:,1),t,y(:,2),t,y(:,3),t,yr(:,1),t,yr(:,2),t,yr(:,3));
legend('$y_1$','$y_2$','$y_3$','$\hat{y}_1$','$\hat{y}_2$','$\hat{y}_3$','interpreter','latex');
xlabel('$t$','interpreter','latex');
ylabel('$y$','interpreter','latex');

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
sigma = deglexord(-1,3,3);
[L,~] = size(sigma);

N = ceil(1.9*L); %N of data points
%N = 50;

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
scatter3(Y(:,1),Y(:,2),Y(:,3),23,'MarkerEdgeColor','g','MarkerFaceColor','y','LineWidth',1.5);

%reconstruct order ideal
eps = 1e-3;
[~, O] = ApproxBM(Y, eps, sigma); %use approximate Buchberger-Moller algorithm

%Use LSM for fitting the equations with the proper coefficients
eta = 1e-2;
H = cell(1,3);
T = cell(1,3);
%reconstruct each equation
for i = 1:3
    V = W(:,i);
    
    [L, ~] = size(O);
    
    E = EvalPoly(eye(L),Y,O);
    h0 = (E'*E)\(E'*V);
    [hi,tau] = delMinorTerms(Y,V,O,eta,h0,1); %get equation and basis
    V0 = EvalPoly(hi,Y,tau); 
    
    H{1,i} = hi;
    T{1,i} = tau;
end

%display equations
prettyABM(H,T)

%simulate results
[~,y] = ode113(@(t,x)oderecon(H,T,t,x),[0:h:Tmax],[0.1,0,-0.1],opts); %solve ODE

figure(1);
plot3(y(:,1),y(:,2),y(:,3),'-');

%compare with analytically derived Lorenz
[~,y] = ode113(@Lorenz2,[0:h:Tmax],[0.1,0,-0.1],opts); %solve ODE

figure(1);
plot3(y(:,1),y(:,2),y(:,3),'-');

legend('full data','data for reconstruction','reconstruction','solution of true equation');

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 