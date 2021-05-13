%% Implementation of Forward
function [alpha, scaleF,loglik] = Forward_DISC(Q, R, pi0, Ob)
% Calculates the Forward Variable
Ns = size(Q,1);
Nt = size(Ob,2);

% step 1: initialization
for j = 1:Ns
    alpha(j,1) = pi0(j)*R(j,Ob(1));
    
end
% normalization
[alpha(:,1), scaleF(1)] = normalise(alpha(:,1));


% step 2: induction (forward)
for t = 2:Nt
    for j = 1:Ns
        sums = 0;
        for i = 1:Ns
            sums = sums+alpha(i,t-1)*Q(i,j);
        end
        alpha(j,t) = sums*R(j,Ob(t));
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
 