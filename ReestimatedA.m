function Abar = ReestimatedA(gamma, xi)
% INPUT
% gamma : smoothing result, P(X_t|Y_{1:T}_
% xi : P(x_t,x_t+1|Y_{1:T})

nTime = size(gamma,2);
nLatentStates = size(gamma,1);

for i = 1:nLatentStates
    denomA = 0;
    for t = 1:nTime-1
        denomA = denomA+gamma(i,t);
    end
    
    for j = 1:nLatentStates
        numA = 0;
        for t = 1:nTime-1
            numA = numA+xi(i,j,t);
        end
%         Abar(i,j) = .001+.999*numA/denomA;
        Abar(i,j) = numA/denomA;
    end
end
