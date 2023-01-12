function [P,A,vmax]=caseTomo_reals_to_P_area(reals,dx,vmax);

if nargin<3, vmax=0.13;end

P=mean(reals<vmax);
Apixel=(reals<vmax).*(dx*dx);
A=sum(Apixel,2);