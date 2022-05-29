function [alw_corr,shd_error,a,bb,aph,ag]=shade_corr_alw(wl,alw_shade,solz)

% correction of the shade error of measured water-leaving albedo using the
% skylight-blocked approach.
% ref: Zhehai Shang, Xiaolong Yu *, Zhongping Lee, A direct measurement system of water-leaving albedo in the field by the skylight-blocked
% approach: Monte Carlo simulations, submitted to Optics Express.

% created by Xiaolong Yu (xlyu@xmu.edu.cn)
 
%% input 
% alw_shade:  measured or simulated water-leaving albedo with shade error 
% wl:         wavelengths match alw_shade, in nm
% solz:       solar zenith angle, in deg
% opt:        input parameters(G-LUT, inherent optical properties of pure
%             seawater)
%% output
 % alw_corr: corrected water-leaving albedo with shade error 
 % a:        total absorption coefficients, in m-1
 % bb:       total backscattering coefficients, in m-1
 % Rrs:      remote sensing reflectance, in sr-1
 
 %% main  

% parameterization of aph
temp=load('a0a1.m');
a0=interp1(temp(:,1),temp(:,2),wl,'PCHIP');
a1=interp1(temp(:,1),temp(:,3),wl,'PCHIP'); 
aw=h2o_iops_Zhh_lee(wl,'a');
bbw=h2o_iops_Zhh_lee(wl,'bb');

%# apply QAA to obtain the initial guess of spectral optimization
% Rrs_temp=alw_shade./pi();
% [a_ini, ~, bbp_ini, eta] = qaa_v6(wl,Rrs_temp,aw,bbw);
% [~,wl550]=min(abs(wl-550));
% [~,wl440]=min(abs(wl-440));
%## initial guess implemented###     
%## five variables: Aph440, Adg440, Bbp555, sdg, eta

bbp550= 0.01; 
P = 0.1;
G = P; 
Sdg = 0.015;
eta=1;

x0=[P,G,bbp550,Sdg,eta];  %initial guess of four variables   
%## spectral optimization ### 
options=optimset('largescale','on','display','iter','tolx',1e-10,'tolfun',1e-12,'Algorithm','active-set','MaxFunEvals',10000,'MaxIter',10000,'Display','off');
LB=[1e-5;1e-5;1e-5;0.008;-0.5];             %lower boundary
UB=[5;5;1;0.02;2.5];                         %upper boundary  
[xx,~]=fmincon(@(x0) shd_error_opt(x0,wl,aw,bbw,solz,a0,a1,alw_shade),x0,[],[],[],[],LB,UB,[],options); 
[alw_corr,~,shd_error,a,bb,aph,ag]=model_alphaw(xx,wl,aw,bbw,solz,a0,a1);
 

end