function [gamma, xi,loglik] = smoother_DISC(pi0, b, R, Q)
% Smoother: P(X_t|Y_{1:T})
% Input
% pi0 : initial probability of X
% b : observation of Y
% GTM : ground truth matrix of Y
% Q : transition marix of X
% variance : variances of multivariate Gaussian

% Output
% gamma : P(x_t|Y_{1:T})
% Xsmoother : chosen X 
% xi : P(x_t,x_{t+1}|Y_{1:T})

[alpha, scaleF,loglik] = Forward_DISC(Q, R, pi0, b);
[beta,gamma] = Backward_DISC(Q, R, b,alpha);

% Baum-Welch's P(x_t,x_{t+1}|Y_{1:T})
xi = ComputeXi_DISC(Q,R,alpha,beta, b);   
xi;

