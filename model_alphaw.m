
function [alphaw,alphaw_shd,shade,a,bb,aph,ag]=model_alphaw(x,wl,aw,bbw,solz,a0,a1)
%  bbp550= bbp_ini(wl550); 
% P = 0.5*(a_ini(wl440)-aw(wl440)); % Aph440
% G = P; % Adg440
% Sdg = 0.015;  
% x0=[P,G,bbp555,Sdg,eta];  %initial guess of four variables  
P=x(1);
G=x(2);
bbp550=x(3);
sdg=x(4);
Y=x(5);

Wn=length(wl);
wl=wl(1:Wn);
aph(1:Wn)=(a0(1:Wn)+a1(1:Wn)*log(P))*P;
ag(1:Wn)=exp(-sdg.*(wl-440))*G;
bbp(1:Wn)=bbp550.*((550./wl).^Y);
a(1:Wn)=aw(1:Wn)+aph(1:Wn)+ag(1:Wn);
bb(1:Wn)=bbw(1:Wn)+bbp(1:Wn);     

[alphaw]=get_alphaw(bbw,bbp,a,solz);
[shade]=get_shderror(a,bb,solz); 
alphaw_shd=alphaw.*(1-shade);
    
end



