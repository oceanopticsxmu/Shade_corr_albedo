function f=optimization_Ew_pre(x,Wn,a,bbw,s0,sza,ntheta,nphi,senz,phi,G_L,Ewsraw)

% Ew=getEw(G_L,bbw,bbp,a,Wn,ntheta,nphi,senz,phi);
Ew=getEw(G_L,bbw,x(1),a,Wn,ntheta,nphi,senz,phi);
shade=getshade(a,(bbw+x(1)),s0,sza,Wn);
Ews=Ew*(1-shade);
temp=abs(Ews-Ewsraw)/Ewsraw;

% [Ew,Ews,shade,a,bb,aph,ag]=opt_cal(x,wl,Wn,aw,bbw,s0,sza,ntheta,nphi,senz,phi,G_L,a0,a1);    
%     temp=geterror(Ews,Ewsraw,Wn);
    f=temp

end