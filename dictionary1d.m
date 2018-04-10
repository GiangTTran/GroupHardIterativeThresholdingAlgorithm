function D = dictionary1(x)
% Dictionary for 1D data
% Copyright: Hayden Schaeffer, Giang Tran, Rachel Ward
% More information can be found at: 
%     H. Schaeffer, G. Tran and R. Ward, "Learning Dynamical Systems and
%     Bifurcation via Group Sparsity", https://arxiv.org/abs/1709.01558

n=length(x);
D = zeros(7,n);

D(1,:) = ones(1,n);
D(2,:) = x;   
D(3,:) = x.^2;   
D(4,:) = x.^3;   
D(5,:) = x.^4;  
D(6,:) = x.^5;  
D(7,:) = x.^6;
