function [Rbar,mu,sigma] = reestimatedRCont(gamma, nLatentStates, nTime, O, M, normY)
%%% R-Matrix
%%% refer: www.ai.mit.edu/~murphyk/Papers/learncg.pdf
gamma2 = reshape(gamma, [nLatentStates 1 nTime]); % gamma2(i,m,t) = gamma(i,t)
m = zeros(O,nLatentStates,M);
op = zeros(O,O,nLatentStates,M);
ip = zeros(nLatentStates,M);
   
    for i=1:nLatentStates
        for k=1:M
            %%% for explanation: do help mixgauss_Mstep
            w = reshape(gamma2(i,k,:), [1 nTime]); % w(t) = w(i,k,t,l)
            wobs = normY .* repmat(w, [O 1]); % wobs(:,t) = w(t) * obs(:,t)
            m(:,i,k) = m(:,i,k) + sum(wobs, 2); % m(:) = sum_t w(t) * obs(:,t)
            op(:,:,i,k) = op(:,:,i,k) + wobs * normY'; % op(:,:) = sum_t w(t) * obs(:,t) * obs(:,t)'
            ip(i,k) = ip(i,k) + sum(sum(wobs .* normY, 2)); % ip = sum_t w(t) * obs(:,t)' * obs(:,t)
        end
    end
    
postmix = sum(gamma,2);
cov_type = 'full';
[mu, sigma] = mixgauss_Mstep(postmix, m, op, ip, 'cov_type', cov_type);
Rbar = mixgauss_prob(normY, mu, sigma);

