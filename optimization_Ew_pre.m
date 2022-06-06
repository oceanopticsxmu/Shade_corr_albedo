% function f=optimization_Ew_pre(x,Wn,a,bbw,s0,sza,ntheta,nphi,senz,phi,G_L,Ewsraw)
function f=optimization_Ew_pre(x,a,bbw,solz,Ewsraw)

% Ew=getEw(G_L,bbw,bbp,a,Wn,ntheta,nphi,senz,phi);
Ew=get_alphaw(bbw,x(1),a,solz);
shade=get_shderror(a,(bbw+x(1)),solz);
Ews=Ew*(1-shade);
temp=abs(Ews-Ewsraw)/Ewsraw;

% [Ew,Ews,shade,a,bb,aph,ag]=opt_cal(x,wl,Wn,aw,bbw,s0,sza,ntheta,nphi,senz,phi,G_L,a0,a1);    
%     temp=geterror(Ews,Ewsraw,Wn);
    f=temp;

end