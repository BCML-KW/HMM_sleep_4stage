function xi = ComputeXi(Q, GTM, alpha, beta, b,variance)
% compute Xi as in paper by Rabiner in Eq. (36)
nLatentStates = size(Q,1);
nTime = size(b,2);
sigma = eye(size(b,1)).*diag(variance);  % unit variance assumption

for t = 1:nTime-1
    for j = 1:nLatentStates
        R(j,t) = mvnpdf(b(:,t)',GTM(j,:),sigma);
    end
    sum = 0;
    for i = 1:nLatentStates
        for j = 1:nLatentStates
            xi(i,j,t) = alpha(i,t)*Q(i,j)*R(j,t)*beta(j,t+1);
            sum = sum+xi(i,j,t);
        end
    end
    
    for i = 1:nLatentStates
        for j = 1:nLatentStates
            xi(i,j,t) = xi(i,j,t)/sum;
        end
    end
end

end