function Qbar = reestimatedQ(pi,E_SiSj,E_Si)
Qbar = [];
currentState = length(pi);
nextState = length(pi);
for i=1:currentState
    for j=1:nextState
        Qbar(i,j)=E_SiSj(i,j)/E_Si(i,1); %%% new
%         Qbar(i,j)=E_SiSj(i,j)/E_Si(1,i); %%% old

    end
end
    