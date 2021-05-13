%% Implementation of Backward
function [beta,gamma] = Backward_DISC(Q, R, Ob,alpha)
% Calculate the Backward Variable
Ns = size(Q,1);
Nt = size(Ob,2);

% step 1: initialization
for j = 1:Ns
    beta(j,Nt) = 1;
end

% gamma (smoother)
gamma(:,Nt) = normalise(alpha(:,Nt).*beta(:,Nt));


% step 2: induction (backward)
for t = Nt-1:-1:1
    for i = 1:Ns
        sums = 0;
        for j = 1:Ns
            sums = sums+Q(i,j)*R(j,Ob(t+1))*beta(j,t+1);
        end
        beta(i,t) = sums;
    end
    % Normalization of beta using scale factor from Forward process
    beta(:,t) = normalise(beta(:,t));
    gamma(:,t) = normalise(alpha(:,t).*beta(:,t));

end

end
