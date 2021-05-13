%% Implementation of Forward
function [alpha, scaleF,loglik] = Forward_CONT(pi0, Q, R)
% Calculates the Forward Variable
[nLatentStates,nTime] = size(R);

% step 1: initialization
% Calculate R matrix, P(Y_t = b(t) | X_t = a(t)), assuming multivariate
% Gaussian distribution likelihood

for j = 1:nLatentStates
    alpha(j,1) = pi0(j)*R(j,1);       % Bayesian
end
% normalization
[alpha(:,1), scaleF(1)] = normalise(alpha(:,1));

% step 2: induction (forward)
for t = 2:nTime
    
    for j = 1:nLatentStates
        sums = 0;
        for i = 1:nLatentStates
            sums = sums+alpha(i,t-1)*Q(i,j);      % One step prediction
        end
        alpha(j,t) = sums*R(j,t);        % Bayesian
    end
    % normalization
    [alpha(:,t), scaleF(t)] = normalise(alpha(:,t));
    
end

if any(scaleF==0)
    loglik = -inf;
else
    loglik = sum(log(scaleF));
end

end
