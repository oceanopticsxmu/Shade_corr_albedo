
function [error]=get_costfunction_error(alphaw_shd,alphaw_shd_raw)
%calculate the error between simulated Rrs and measured Rrs     
%error = {mean[(raw-sim)^2)]}^0.5/mean[(raw)]
Wn=length(alphaw_shd);
temp1=sum(alphaw_shd_raw(1:Wn));  
temp2=sum((alphaw_shd-alphaw_shd_raw).^2);
temp3=temp1/Wn;
temp4=temp2/Wn;
%error=((temp1+temp3)^0.5)/(temp2+temp4);
error=(temp4^0.5)/temp3;
    
end