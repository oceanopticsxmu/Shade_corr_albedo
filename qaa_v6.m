function [a, bb, bbp, Y] = qaa_v6(wl,Rrs,aw,bbw)
% updated Qusi-analytic algorithm (QAA-v6)
% Function call: bbw.m (backscattering coefficient of pure sea water)
% Inputs:
%   Rrs: remote-sensing reflectance
%   wl: wavelength
% Outputs:
%   a, bb: total absorbtion and backscattering
%   anw: nonwater absorption, anw= a - aw
%   bbp: particulate backscattering
%   aph: phytoplankton abs
%   acdm: CDM absorption
% Functions to call
%   h2o_iops_Zhh_lee
% *************************************************************
% Citation: Lee Z. P., IOCCG, online link, quasi-analytical algorithm (QAA-v6); not published
% modified Feb 3, 2021, by Xiaolong Yu, using pure seawater absorption (Lee et al., 2005) and
% backscattering (Zhang et al., 2009)
% ***************************************************************************
% aw bbw
nwl=length(wl);


% id412 = find(abs(wl-412)==min(abs(wl-412))); 
id443 = find(abs(wl-443)==min(abs(wl-443))); 
id490 = find(abs(wl-490)==min(abs(wl-490))); 
id555 = find(abs(wl-555)==min(abs(wl-555))); 
id670 = find(abs(wl-670)==min(abs(wl-670)));  
id555=id555(1);

%%%% Step 1
rrs(1:nwl) = Rrs./(0.52+1.7*Rrs);

g0 = 0.089;
g1 = 0.125;
u  = (-g0 + (g0^2 + 4*g1*rrs).^0.5)/(2*g1);


%%%% Step 2

if Rrs(id670) < 0.0015     
    
    wl_ref = wl(id555);
    id_ref = id555;     
    %%%%
%     rrs_ref = rrs(id_ref);   
    ki = log10((rrs(id443)+rrs(id490))/(rrs(id_ref) + 5*rrs(id670)*rrs(id670)/rrs(id490)));
%     a_ref = aw_lee2015(wl_ref) + 10.^(-1.146 - 1.366 * ki - 0.469*ki.^2);
    a_ref = aw(id_ref) + 10.^(-1.146 - 1.366 * ki - 0.469*ki.^2);
else
    wl_ref = 670;
    id_ref = id670;
    a_ref = aw(id_ref) + 0.39 * (rrs(id670)/(rrs(id443) + rrs(id490)))^1.14;
%     a_ref = aw_lee2015(wl_ref) + 0.39 * (rrs(id670)/(rrs(id443) + rrs(id490)))^1.14;
end


%%%% Step 3

% bbp_ref = u(id_ref).*a_ref/(1-u(id_ref)) - 0.5 * h2o_iops(wl_ref,'b');
bbp_ref = u(id_ref).*a_ref/(1-u(id_ref)) - bbw(id_ref);
%%%% Step 4
Y = 2.0*(1-1.2*exp(-0.9*rrs(id443)/rrs(id555))); 

%%%% Step 5 & Step 6
% for i = 1 : length(wl)
%    bbp(i,1) = bbp_ref *(wl_ref/wl(i))^Y;
%    a(i,1)   = (1-u(i)).*(bbw(i) + bbp(i))./u(i);    
% %    a(i,1)   = (1-u(i)).*(0.5*h2o_iops(wl(i),'b') + bbp(i))./u(i);  
% end

bbp(1:nwl)=bbp_ref *(wl_ref./wl).^Y;
bb(1:nwl) = bbp + bbw; 
a(1:nwl) = (1-u).*bb./u;   

%% decomposition to component absorptions
% %%%% Step 7  (an empirical relationship)
% zeta = 0.74 + 0.2 / (0.8 + rrs(id443)/rrs(id555));  % zeta = a_ph(410)/a_ph(440);
% 
% %%%% Step 8
% S  = 0.015 + 0.002/(0.6+rrs(id443)/rrs(id555));
% xi = exp(S*(443-412));                              % xi = ag(410)/ag(440);
% 
% %%%% Step 9
% acdm443 = (a(id412)-zeta*a(id443))/(xi-zeta) - (aw_lee2015(412)-zeta*aw_lee2015(443))/(xi-zeta); 
% 
% acdm = acdm443 * exp(-S*(wl - 443));
% 
% 
% %%%% Step 10
% aph = a' - acdm - aw_lee2015(wl);
% 
% anw = a' - aw_lee2015(wl); 

end


