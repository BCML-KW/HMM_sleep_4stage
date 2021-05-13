function [pi0,Q,Rupdate] = EM_CONT(pi0,b,Rall,R,Q,iterN)
% INPUT
% pi0 : intial probability of X
% b : observation
% R : P(Y|X)
% Q : transition of X, P(X_t|X_{t+1})
% V : variables of Y
% iterN : iteration number

[nLatentStates,nTime,M] = size(Rall);
[O,~] = size(b);

% initial parameters for converging test
previous_loglik = -inf;
threshold = 1e-4;   

converged = 0;
num_iter = 1;

while (num_iter <= iterN) & ~converged    
    %% pi and Q matrix update
    % E Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Discrete for Q (transition) matrix
    time = 1:nTime;    % manipulate discrete function
    [gamma, xi,current_loglik] = smoother_DISC(pi0, time, R, Q);
    disp(['iteration = ',num2str(num_iter), '    log(P(Y|parameter) = ',num2str(current_loglik)]);
    
    [converged, decrease] = em_converged(current_loglik, previous_loglik, threshold);
    previous_loglik = current_loglik;
    
    %%%%%% PLOT Convergence of Log-Likelihood %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(num_iter,current_loglik,'b*');
    title('Convergence of Log-Likelihood');xlabel('Number of Iterations');ylabel('Log-Likehood Value');
    hold on; 
    
    % M Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % update transition matrix
    newQ = ReestimatedA(gamma, xi); % Updated transition matrix
    Q = newQ;
    
    % update pi matrix
    pi0 = gamma(:,1)';
    
    %% R matrix update    
    % See www.ai.mit.edu/~murphyk/Papers/learncg.pdf.
    gamma2 = reshape(gamma, [nLatentStates 1 nTime]); % gamma2(i,m,t) = gamma(i,t)
    
    m = zeros(O,nLatentStates,M);
    op = zeros(O,O,nLatentStates,M);
    ip = zeros(nLatentStates,M);
    
    for i=1:nLatentStates
        for k=1:M
            w = reshape(gamma2(i,k,:), [1 nTime]); % w(t) = w(i,k,t,l)
            wobs = b .* repmat(w, [O 1]); % wobs(:,t) = w(t) * obs(:,t)
            m(:,i,k) = m(:,i,k) + sum(wobs, 2); % m(:) = sum_t w(t) * obs(:,t)
            op(:,:,i,k) = op(:,:,i,k) + wobs * b'; % op(:,:) = sum_t w(t) * obs(:,t) * obs(:,t)'
            ip(i,k) = ip(i,k) + sum(sum(wobs .* b, 2)); % ip = sum_t w(t) * obs(:,t)' * obs(:,t)
        end
    end
    
    postmix = sum(gamma,2);
    cov_type = 'full';
    [mu, sigma] = mixgauss_Mstep(postmix, m, op, ip, 'cov_type', cov_type);
        
%     R = mixgauss_prob(b, mu, sigma);
    R = R_matrix(b,permute(mu,[2 1]),permute(sigma,[3 1 2]));
    num_iter = num_iter+1;
end


Rupdate.mu = permute(mu,[2 1 3]);
Rupdate.sigma = squeeze(permute(sigma,[3 1 2 4]));
Rupdate.pi0 = pi0;