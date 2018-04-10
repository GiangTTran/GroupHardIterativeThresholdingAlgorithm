% Group Hard Iterative Thresholding Algorithm on Lorenz system. Test 5 data
% sets.

% Note: The $\ell_0$-model is the group hard-iterative thresholding algorithm
% applied to only one data set at a time.

% Copyright: Hayden Schaeffer and Giang Tran and Rachel Ward 2017
% More information can be found at: 
%     H. Schaeffer, G. Tran and R. Ward, "Learning Dynamical Systems and
%     Bifurcation via Group Sparsity", https://arxiv.org/abs/1709.01558
close all; clear all; clc

%% Parameters
T = 20.0;
dt = 0.005;
sigma_noise = dt;% noise level in data
m = 5; % number of data sets

sigma = 10.0; % a coefficient in the Lorenz system
beta = 8.0/3.0; % a coefficient in the Lorenz system

c1 = -1.0; % original Lorenz system - chaos
U01 = [-8.0 7.0 27.0]; % initial condition for parameter c1


c2 = 4.70; % original Lorenz system - chaos
U02 = [0.0 -0.01 9.0]; % initial condition for parameter c1

c3 = 6.9; % 
U03 = [1.0 2.0 1.0];

c4 = 7.075; % limit cycle - 
U04 = [1.0 1.0 2.0]; % initial condition for parameter c3

c5 = 7.730;
U05 = [2.0 1.0 -5.0];

thres = 2.1; % threshold parameter of the group hard-iterative thresholding algorithm
numIter = 10; % number of iterations of the group hard-iterative thresholding algorithm

%% Generate Data by solving the ODE system with ode45
%Input 1
c = c1;
U0 = U01;

myfunc = @(t,x) [sigma*(x(2)-x(1))   ;  - x(1).*x(3)+(24-4*c)*x(1) + c*x(2) ;   x(1).*x(2)-beta*x(3)];
[tvec X1] = ode45(myfunc,[0:dt:T],U0,odeset('RelTol',1e-12, 'AbsTol',[1e-12 1e-12 1e-12]));

%Input 2
c = c2;
U0 = U02;

myfunc = @(t,x) [sigma*(x(2)-x(1))   ;  - x(1).*x(3)+(24-4*c)*x(1) + c*x(2) ;   x(1).*x(2)-beta*x(3)];
[tvec X2] = ode45(myfunc,[0:dt:T],U0,odeset('RelTol',1e-12, 'AbsTol',[1e-12 1e-12 1e-12]));


%Input 3
c = c3;
U0 = U03;

myfunc = @(t,x) [sigma*(x(2)-x(1))   ;  - x(1).*x(3)+(24-4*c)*x(1) + c*x(2) ;   x(1).*x(2)-beta*x(3)];
[tvec X3] = ode45(myfunc,[0:dt:T],U0,odeset('RelTol',1e-12, 'AbsTol',[1e-12 1e-12 1e-12]));

%Input 4
c = c4;
U0 = U04;

myfunc = @(t,x) [sigma*(x(2)-x(1))   ;  - x(1).*x(3)+(24-4*c)*x(1) + c*x(2) ;   x(1).*x(2)-beta*x(3)];
[tvec X4] = ode45(myfunc,[0:dt:T],U0,odeset('RelTol',1e-12, 'AbsTol',[1e-12 1e-12 1e-12]));


%Input 5
c = c5;
U0 = U05;

myfunc = @(t,x) [sigma*(x(2)-x(1))   ;  - x(1).*x(3)+(24-4*c)*x(1) + c*x(2) ;   x(1).*x(2)-beta*x(3)];
[tvec X5] = ode45(myfunc,[0:dt:T],U0,odeset('RelTol',1e-12, 'AbsTol',[1e-12 1e-12 1e-12]));

% Generate noisy data
X1noise = X1 + sigma_noise*max(abs(X1(:)))*randn(size(X1));
X2noise = X2 + sigma_noise*max(abs(X2(:)))*randn(size(X2));
X3noise = X3 + sigma_noise*max(abs(X3(:)))*randn(size(X3));
X4noise = X4 + sigma_noise*max(abs(X4(:)))*randn(size(X4));
X5noise = X5 + sigma_noise*max(abs(X5(:)))*randn(size(X5));

% Approximate velocity
v1 = time_derivative(X1noise,dt,2); 
v2 = time_derivative(X2noise,dt,2); 
v3 = time_derivative(X3noise,dt,2); 
v4 = time_derivative(X4noise,dt,2); 
v5 = time_derivative(X5noise,dt,2); 

% Compute dictionary matrices 
D1 = dictionary3(X1noise(2:end-1,:)); % Dictionary - v1 and D1 have the same #columns

D2 = dictionary3(X2noise(2:end-1,:)); % Dictionary - v1 and D1 have the same #columns

D3 = dictionary3(X3noise(2:end-1,:)); % Dictionary - v1 and D1 have the same #columns

D4 = dictionary3(X4noise(2:end-1,:)); % Dictionary - v1 and D1 have the same #columns

D5 = dictionary3(X5noise(2:end-1,:)); % Dictionary - v1 and D1 have the same #columns

%% Visualization of input data
figure('name','Noisy Data Sets');
subplot(231); plot3(X1noise(:,1),X1noise(:,2),X1noise(:,3),'LineWidth',1.5);axis image; view(27,16); 
subplot(232); plot3(X2noise(:,1),X2noise(:,2),X2noise(:,3),'LineWidth',1.5);axis image; view(27,16); 
subplot(233); plot3(X3noise(:,1),X3noise(:,2),X3noise(:,3),'LineWidth',1.5);axis image; view(27,16); 
subplot(234); plot3(X4noise(:,1),X4noise(:,2),X4noise(:,3),'LineWidth',1.5);axis image; view(27,16); 
subplot(235); plot3(X5noise(:,1),X5noise(:,2),X5noise(:,3),'LineWidth',1.5);axis image; view(27,16); 

figure('name','Velocity Sets');
subplot(231); plot3(v1(:,1),v1(:,2),v1(:,3),'.');axis image; view(27,16); 
subplot(232); plot3(v2(:,1),v2(:,2),v2(:,3),'.');axis image; view(27,16); 
subplot(233); plot3(v3(:,1),v3(:,2),v3(:,3),'.');axis image; view(27,16); 
subplot(234); plot3(v4(:,1),v4(:,2),v4(:,3),'.');axis image; view(27,16); 
subplot(235); plot3(v5(:,1),v5(:,2),v5(:,3),'.');axis image; view(27,16); 

%% Parameters + Input for the group hard-iterative thresholding algorithms
% Ground truth coefficients of dy/dt = - xz+(24-4*c)*x + c*y
indexTrue = [2,3,7];
C_true = zeros(size(D1,2),m);
C_true(indexTrue,1) = [24-4*c1,c1,-1]; % coefficients of dy/dt from Data 1
C_true(indexTrue,2) = [24-4*c2,c2,-1]; % coefficients of dy/dt from Data 2
C_true(indexTrue,3) = [24-4*c3,c3,-1]; % coefficients of dy/dt from Data 3
C_true(indexTrue,4) = [24-4*c4,c4,-1]; % coefficients of dy/dt from Data 4
C_true(indexTrue,5) = [24-4*c5,c5,-1]; % coefficients of dy/dt from Data 5

% time step
ds = 1/max([norm(D1'*D1,'fro'),norm(D2'*D2,'fro'),norm(D3'*D3,'fro'), norm(D4'*D4,'fro'), norm(D5'*D5,'fro')]);

% 
Dgroup = cell(m,1);
Dgroup{1} = D1;
Dgroup{2} = D2;
Dgroup{3} = D3;
Dgroup{4} = D4;
Dgroup{5} = D5;

% recover dx2/dt
rhsgroup = cell(m,1);
rhsgroup{1} = v1(:,2); 
rhsgroup{2} = v2(:,2); 
rhsgroup{3} = v3(:,2); 
rhsgroup{4} = v4(:,2); 
rhsgroup{5} = v5(:,2); 

% Initialization
Cinitial = ones(size(D1,2),m);

%% Group Hard-Iterative Thresholding Algorithm
C = group_hard_iterative_thresholding(Dgroup,rhsgroup,ds,Cinitial,thres,numIter);

%% Compare with the ground truth

% Check whether C_true and C have the same support set
if (length(find((abs(C)>0) - (abs(C_true)>0)))==0)
    fprintf('The group hard-iterative algorithm correctly identifies all terms in the governing equation')
else
    fprintf('The group hard-iterative algorithm doesnot identify correctly all terms in the governing equation')
end

% Check relative error
error_rel = max(abs(C-C_true),[],1)./max(abs(C_true),[],1);
fprintf('The relative errors in the coefficients of the govering equation associated with each data set are\n');
error_rel

return
