function [ ] = ChangeVar( Var, var )
%This function changes the value of the variable of interest
global S1 S2 S3
    if strcmp(Var,'GM-CSF')
        S1=var;
    elseif strcmp(Var,'M-CSF')
        S2=var;
    elseif strcmp(Var,'G-CSF')
        S3=var;
    else
        error('incorrect Var')
    end
end

