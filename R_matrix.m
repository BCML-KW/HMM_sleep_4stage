function R = R_matrix(b,mu,sigma)

[nLatentStates,nO] = size(mu);
[~, nTime] = size(b);
%% HMM Toolbox
R = mixgauss_prob(b, permute(mu,[2 1]), squeeze(permute(sigma,[2 3 1])));

%% Matlab code
% for t = 1:nTime
%     for j = 1:nLatentStates
%         R(j,t) = mvnpdf(b(:,t)',squeeze(mu(j,:)),squeeze(sigma(j,:,:)));
%     end
% end