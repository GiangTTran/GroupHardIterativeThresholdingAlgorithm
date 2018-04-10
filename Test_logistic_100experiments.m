% Group-Sparse Regression on 1d ODE with multi-data
%     xdot = alpha*x*(1-x), 
% where alpha is the bifurcation parameter.
% Output: the probability of recovering correct terms and the relative error ranges

% Copyright: Hayden Schaeffer, Giang Tran, Rachel Ward
% More information can be found at: 
%     H. Schaeffer, G. Tran and R. Ward, "Learning Dynamical Systems and
%     Bifurcation via Group Sparsity", https://arxiv.org/abs/1709.01558

% Add noise to data and estimate derivative from noisy data

close all; clear all; clc
p = 6; % maximal degree of monomials
T = 50.0;
dt = 0.005;

sigma_noise = 2*dt;
numIter = 10; % number of iterations of the group hard-iterative thresholding algorithm
thres = 2e-3;
thres1 = 5e-6;
thres2 = 2e-3;

m = 2; % number of data sets
maxIter = 100; % number of experiments

count = 0; % probability of recover the governing equation for both data sets using group sparsity
count11 = 0; % probability of recover the governing equation for data set 1 using l_0-model
count12 = 0; % probability of recover the governing equation for data set 2 using l_0-model

error_rel = zeros(maxIter,2); % relative error in the recovered coefficient vector for both data sets using group sparsity

%Input 1
x0 = 0.01; 
rhs = @(t,x) (0.05*x.*(1-x));
[~, x1_clean] = ode45(rhs,0:dt:T,x0,odeset('RelTol',1e-12,'AbsTol',1e-12));

%Input 2
x0 = 0.01; 
rhs = @(t,x) (0.23*x.*(1-x));
[tvec, x2_clean] = ode45(rhs,0:dt:T,x0,odeset('RelTol',1e-12,'AbsTol',1e-12));

indexTrue =[2,3];
nbar = p+1; % number of basis terms
C_true = zeros(nbar,m);
C_true(indexTrue,1) = [0.05,-0.05]; % coefficients of dy/dt from Data 1
C_true(indexTrue,2) = [0.23,-0.23]; % coefficients of dy/dt from Data 2


for iter = 1:maxIter
    % Generate noise data
    x1 = x1_clean + sigma_noise*max(abs(x1_clean))*randn(size(x1_clean));
    x2 = x2_clean +  sigma_noise*max(abs(x2_clean))*randn(size(x2_clean));

    % Approximate velocity    
    v1 = time_derivative1(x1,dt,2);
    v2 = time_derivative1(x2,dt,2);
    
    % Compute dictionaries
    D1 = dictionary1d(x1)';
    D2 = dictionary1d(x2)';

    Dgroup = cell(m,1);
    Dgroup{1} = D1;
    Dgroup{2} = D2;
    
    rhsgroup = cell(m,1);
    rhsgroup{1} = v1; 
    rhsgroup{2} = v2;
    
     % Initialization
    Cinitial = zeros(nbar,m);
    
    ds = 1/max([norm(D1'*D1,'fro'),norm(D2'*D2,'fro')]);
    %% Group Hard-Iterative Thresholding Algorithm
    C = group_hard_iterative_thresholding(Dgroup,rhsgroup,ds,Cinitial,thres,numIter);
    %% Compare with the ground truth

    % Check whether C_true and C have the same support set
    if (length(find((abs(C)>0) - (abs(C_true)>0)))==0)
        count = count+1;
        error_rel(iter,:) = max(abs(C-C_true),[],1)./max(abs(C_true),[],1);
    end

    %% Main -- Only Data 1 -- Hard-Iterative Thresholding Solution
    %Parameter for sparsity
    C11 = zeros(size(D1,2),1); %Initialize

    %Step size
    A1 = D1'*D1;
    ds1 = 1/max(norm(A1,'fro'));
    b1 = D1'*v1;
    
    % Hard Iterative Thresholding
    for j = 1: 10
        C11 = C11-ds1*(A1*C11-b1);
        
        C = C11;
        C_sum = sqrt(sum( C.^2,2));

        Index = C_sum < thres1;
        NIndex = C_sum > thres1;
        
        C11(Index) = 0;
        C11(NIndex) =    D1(:,NIndex)\v1(:);
    end

     if (isempty(setdiff(indexTrue,find(C11))) && isempty(setdiff(find(C11),indexTrue)))
            count11 = count11+1;
     end


    %% Main -- Only Data 2 -- Hard-Iterative Thresholding Solution
    C12 = zeros(size(D2,2),1);

    %Step size
    A2 = D2'*D2;
    ds2 = 1/max(norm(A2,'fro'));
    b2 = D2'*v2;
    
    for j = 1: 10
        C12 = C12-ds2*(A2*C12-b2);
        C = C12;

        C_sum = sqrt(sum( C.^2,2));
        Index = C_sum < thres2;
        NIndex = C_sum>thres2;

        C12(Index) = 0;
        C12(NIndex) =    D2(:,NIndex)\v2(:);

    end
    if (isempty(setdiff(indexTrue,find(C12))) && isempty(setdiff(find(C12),indexTrue)))
            count12 = count12+1;
    end

end
% Probability of recovering correct terms using group sparsity
fprintf(['Probability of recovering correct terms (out of ',num2str(maxIter), ') using l^{2,0} is ',num2str(count),'\n'])

% % Relative error ranges
% tmp = error_rel;
% tmp(tmp==0) = NaN;
% error_min = min(tmp,[],1);
% error_max = max(error_rel,[],1);
% fprintf('Relative error ranges are \n ')
% [error_min;error_max]

fprintf(['Probability of recovering correct terms (out of ',num2str(maxIter), ') for data set 1 and data set 2 using l^0 are ',num2str(count11),' and ', num2str(count12),'\n']); 

% Visualization of an experiment
figure; 
subplot(221); plot(tvec,x1,'LineWidth',1.5); set(gca,'FontSize',18);axis tight


subplot(223); plot(tvec,v1,'.'); set(gca,'FontSize',18); axis tight

subplot(222);
plot(tvec,x2,'LineWidth',1.5); set(gca,'FontSize',18); axis tight

subplot(224);
plot(tvec,v2,'.'); set(gca,'FontSize',18); axis tight

