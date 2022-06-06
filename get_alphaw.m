function [alphaw]=get_alphaw(bbw,bbp,a,solz)
%calculate alpha_w

senz=[0,10,20,30,40,50,60,70,80,87.5];
phi=[0,15,30,45,60,75,90,105,120,135,150,165,180]; 
ntheta=length(senz);
nphi=length(phi);
 
%##get alpha_w ##################
Gfilename='G_LUT.mat';
load(Gfilename);  
G0=get_G0(solz,G_LUT,senz,phi);        
k=a+bbw+bbp;
Wn=length(bbw);

for i = 1: ntheta
    for j=1:nphi
        Angular_Rrs(i,j,1:Wn)=(G0(i,j,1)+G0(i,j,2)*bbw./k).*bbw./k+(G0(i,j,3)+G0(i,j,4)*bbp./k).*bbp./k; 
    end
end
         
for k = 1:length(bbw)
    A_Rrs=Angular_Rrs(:,:,k);
    alphaw(k)=getEw_interp2(A_Rrs,senz,phi);
end
     
end
