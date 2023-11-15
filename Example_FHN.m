%Reconstruction of FitzHugh-Nagumo system using ABM algorithm and delMinorTerms
close all

Tmax = 40; %overall time
h = 0.02;

Ttr = 0; %transient
[t,y] = ode45(@FHN,[0:h:Tmax],[1,0]); %solve ODE

eps = 1e-2; %tolerance for ABM algorithm
eta = 1e-2; %tolerance for LSM algorithm

derivativecalc = 1; % Set 1 to find derivatives analytically, set 0 to find derivatives numerically. Check out how this affects accuracy

if derivativecalc
    w = transpose(FHN(0,transpose(y)));
else
    w = diff4(y,t); 
end

y = y(1:end-1,:);

figure(1); 
plot(y(:,1),y(:,2));
xlabel('\itx');
ylabel('\ity');

%get uniformly distributed points
N = 20; %data points
M = 2; %dim
[Ns, ~] = size(y);
W = zeros(N,M);
Y = zeros(N,M);

Ntr = Ttr/h; %number of points in transient

for i = 1:N
    id = ceil(rand*(Ns - Ntr) + Ntr); %number of data point
    W(i,:) = w(id,:); %derivatives
    Y(i,:) = y(id,:); %Y
end

%plot points
figure(1); hold on
scatter(Y(:,1),Y(:,2),23,'MarkerEdgeColor','g','MarkerFaceColor','y','LineWidth',1.5);

[H,T] = PolyRegression(Y,W,0,3,eta,eps); %reconstruct H and T with one function

%simulate results
[~,y] = ode45(@(t,x)oderecon(H,T,t,x),0:h:Tmax,[1,0]); %solve ODE
figure(1);
plot(y(:,1),y(:,2));

%display equations
prettyABM(H,T)

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 