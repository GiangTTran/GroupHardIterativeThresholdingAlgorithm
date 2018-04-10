function C = group_hard_iterative_thresholding(Dgroup,rhsgroup,ds,Cinitial,thres,numIter)
% Solve 
%
% min_{c1,...,cm} \|D1 * c1 - b1\|_2^2 + ... \|Dm * cm - bm\|_2^2 + gamma * \|[c1,..., cm]\|_{2,0}  
%
% where \|A\|_{2,0} = number of non-zero rows of the matrix A.

% Input:
%      Dgroup = {D1,...,Dm} where Dk is of size l_k x nbar
%      rhsgroup = {v1,...,vm} where vk is of size l_k x 1
%      Cinitial = Initial values of the unknown matrix C=[c1... cm]
%      thres = thresholding parameter (=\sqrt(gamma))
%      numIter = number of group hard-iterative algorithm
% Output:
%      C = matrix of size nbar x m
%
% Copyright: Hayden Schaeffer, Giang Tran, Rachel Ward
% More information can be found at: 
%     H. Schaeffer, G. Tran and R. Ward, "Learning Dynamical Systems and
%     Bifurcation via Group Sparsity", https://arxiv.org/abs/1709.01558

    m = size(Dgroup,1);
    C = Cinitial;
        
    for j = 1:numIter
        % gradient-descent step
        for i = 1:m
            ctemp = C(:,i);
            Dtemp = Dgroup{i};
            A = (Dtemp)'*Dtemp;
            b = Dtemp'*rhsgroup{i};
            ctemp = ctemp - ds*(A*ctemp - b);
            Cthres(:,i) = ctemp(:);
        end

        C_sum = sqrt(sum( Cthres.^2,2));

        Index = C_sum < (thres);
        NIndex = C_sum>= (thres);

        for i = 1:m
            C(Index,i) = 0;
            Dtemp = Dgroup{i};
            C(NIndex,i) = Dtemp(:,NIndex)\rhsgroup{i};
        end
    end
return;
