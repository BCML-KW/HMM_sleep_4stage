function Bbar = ReestimatedB(gamma, b, V)
% INPUT
% gamma : smoothing result, P(X_t|Y_{1:T})
% b : Observation matrix
% V : variables of Y

nTime = size(gamma,2);
nLatentStates = size(gamma,1);

% Calculate the reestimated B
for i = 1:nLatentStates
    denomA = 0;
    for t = 1:nTime-1
        % Expected number of times in state i
        denomA = denomA+gamma(i,t);
    end
    
    denomB = denomA+gamma(i,nTime);
    for k = 1:length(V)
        numB = 0;
        for t = 1:nTime
            if b(t) == k
                % Expected number of times in state i and observing symbol V(k) 
                numB = numB + gamma(i,t);
            end
        end
        Bbar(i,k) = .001+.999*numB/denomB;
%         Bbar(i,k) = numB/denomB;
    end
end



