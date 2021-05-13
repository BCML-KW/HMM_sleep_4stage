%% Implementaion of Gamma
function gamma_CONT = ComputeGamma_CONT(alpha, beta, Rall,R)
% INPUT
% alpha : forward
% beta : backward
% R : pdf of Gaussian mixture
% OUTPUT
% gamma_CONT : # of states X # of time samples X # of Gaussian mixtures

% computes Gamma in paper by Rabiner, P(X_t|Y_{1:T}),parameter lambda in Rabiners' paper)
[nLatentStates,nTime,M] = size(Rall);

% Discrete Gamma
for t = 1:nTime
    denom = 0;
    for j = 1:nLatentStates
        gamma(j,t) = alpha(j,t)*beta(j,t);      %(Eq. 27)
        denom = denom + gamma(j,t);
    end
    for j = 1:nLatentStates
        gamma(j,t) = gamma(j,t)/(denom);
    end
end
clear denom;
    

for t = 1:nTime
    for j = 1:nLatentStates
        for k = 1:M
            gamma_CONT(j,t,k) = gamma(j,t)*squeeze(Rall(j,t,k))/R(j,t);    
        end
    end       
end

end
