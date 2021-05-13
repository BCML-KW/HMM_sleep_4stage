function smoother_CONT(pi0, b, Q, R,c)
% Smoother: P(X_t|Y_{1:T})
% Input
% pi0 : initial probability of X
% b : observation of Y
% Q : transition marix of X
% R : P(Y|X)

% Output


[alpha, scaleF] = Forward_CONT(pi0, b, Q, R);
beta = Backward_CONT(Q, R, scaleF);
gamma = ComputeGamma_CONT(alpha,beta,R,c);   

% P(Y|X)

