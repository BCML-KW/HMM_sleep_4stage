%% Implementaion of xi
function xi = ComputeXi_DISC(Q, R, alpha, beta, Ob)
% compute Xi as in paper by Rabiner in Eq. (36)
Ns = size(Q,1);
Nt = size(Ob,2);

for t = 1:Nt-1
    sum = 0;
    for i = 1:Ns
        for j = 1:Ns
            xi(i,j,t) = alpha(i,t)*Q(i,j)*R(j,Ob(t+1))*beta(j,t+1);
            sum = sum+xi(i,j,t);
        end
    end
    
    for i = 1:Ns
        for j = 1:Ns
            xi(i,j,t) = xi(i,j,t)/sum;
        end
    end
end

end