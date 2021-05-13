function [gamma, XImatrix, loglik, E_Si, E_SiSj] = smoother_Cont_Sleep(pi,Q,R,Y)
% Smoother: P(X_t|Y_{1:T},Theta)
% Input
% pi = Initial Probability of X
% Q  = Transition Matrix of X
% R = P(Y|X)
% Y = observation values matrix


% Output
% gamma : P(X_t|Y_{1:T})
% XImatrix : P(X_t,X_{t+1}|Y_{1:T})
% E_Si : Expected # transitions from current state Si
%        = sum(gamma_t{i}) from t=1:T-1
% E_SiSj : Expected # transitions from current state Si to next state Sj
%          = sum(XImatrix_t{i,j}) from t=1:T-1

T = length(Y);
nStates = length(pi);

[alpha,loglik] = Fwd_Cont_Sleep(pi,Q,R,Y);
[beta] = Bkwd_Cont_Sleep(pi,Q,R,Y);

% [~,indxc] = find(sum(alpha,1)==0);
% alpha(:,indxc) = alpha(:,indxc)+10^-200*ones(size(alpha(:,indxc)));
% 
% 
% [~,indxcb] = find(sum(beta,1)==0);
% beta(:,indxcb) = beta(:,indxcb)+10^-200*ones(size(beta(:,indxcb)));


[gamma, E_Si] = Gamma_Disc_Sleep(alpha,beta,T,nStates); 
[XImatrix,E_SiSj] = XI_Cont_Sleep(alpha,beta,Q,R,T,nStates);
