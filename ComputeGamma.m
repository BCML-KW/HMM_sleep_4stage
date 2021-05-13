%% Implementaion of Gamma
function [gamma Xpred] = ComputeGamma(alpha, beta)
% computes Gamma in paper by Rabiner, P(X_t|Y_{1:T}),parameter lambda in Rabiners' paper)
nLatentStates = size(alpha,1);
nTime = size(alpha,2);

for t = 1:nTime
    denom = 0;
    for j = 1:nLatentStates
        gamma(j,t) = alpha(j,t)*beta(j,t);      %(Eq. 27)
        denom = denom + gamma(j,t);
    end
    for i = 1:nLatentStates
        gamma(i,t) = gamma(i,t)/(denom);
    end
end

for i = 1:nTime
    [tmp Xpred(i)] = max(gamma(:,i));
end
end
