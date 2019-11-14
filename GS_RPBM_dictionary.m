function [sig,Z,T]=GS_RPBM_dictionary(I)
t     = I.t';
dsize = I.dictionarysize;
dnum  = I.dictionarynumber;

%% Generate Dictionary
for ic=1:dnum
Z(:,ic)  = rand(dsize,1).*5;
T(:,ic)  = rand(dsize,1).*1000;
end

parfor ic=1:dnum
    ic
    for ib=1:dsize
        sig(:,1,ib,ic)=get_Dt_2015(t./T(ib,ic),Z(ib,ic));
    end
end

end

