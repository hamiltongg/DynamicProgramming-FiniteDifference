function F=c_op(c1,a,g1,g2,Y)
%{ 
The goal is to find "c1" which comes from nonlinear function
Date: July 2023
Author: Hamilton Galindo Gil
%}

%% Obtaining c1
%{
c1 = (a^(-1/g1))*(Y - c1)^(g2/g1)
%}

F = c1 - (a^(-1/g1))*(Y - c1)^(g2/g1);
