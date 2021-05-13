function X_map = MAPSEQ_Estimation(b,a,pi0,Q,R)
% INPUT
% b : Observation matrix (feature # X sample #)
% a : latent (hidden) variable X
% B10 : initial probability pi
% Q : transition matrix
% GTM : ground truth matrix for the features (latent variable # X feature #)

%% Maximum a Posteriori (MAP) Sequence Estimation
nLatentStates = length(a);
for t = 1:size(b,2)+1
    if t == 1   % i = 1
        for j = 1:length(a)
            for jj = 1:length(a)
                % d_1(a_0,a_1) in the lecture note, d(t,j,jj) - t:time index, j: state at t-1, jj:state at t
                d(t,j,jj) = -1*log(pi0(jj))-log(R(jj,t));
            end
        end
    elseif t == length(b)+1 % i = n+1
        d(t,:,:) = 0;
    else
        for j = 1:length(a)
            for jj = 1:length(a)
                % d_i(a_{t-1},a_t) in the lecture note, d(t,j,jj) - t:time index, j: state at t-1, jj:state at t
                d(t,j,jj) = -1*log(Q(j,jj))-log(R(jj,t));
            end
        end
    end
end

for t = size(b,2)+1:-1:1
    if t == size(b,2)+1
        for j = 1:length(a)
            J(t, j) = 0;        % J(t,j) - t:time index, j:state index
        end
    else
        for j = 1:length(a)
            [J(t, j) a_map(t,j)] = min(squeeze(d(t,j,:))+J(t+1,:)');    % a_map is a^*_n(a_{n-1}) in the lecture note
        end
    end
end

for t = 1:size(b,2)  
    if t == 1
        X_map(t) = a_map(1,1);
    else
        X_map(t) = a_map(t,X_map(t-1));
    end
end

