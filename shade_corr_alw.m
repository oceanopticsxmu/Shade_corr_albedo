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

%% Apply QAA to help the initial guess and opt limit
%% the index of the wavelength need to be modified if spectral information changes
    s0(1:6)=[0.0137,0.0759,0.1328,0.2298,5,0.6874];
    senz=[0,10,20,30,40,50,60,70,80,87.5];
    phi=[0,15,30,45,60,75,90,105,120,135,150,165,180]; 
    ntheta=length(senz);
    nphi=length(phi);

    id_ref=4; %reference 550nm
    alphaw_shade=alw_shade;
    [~,Wn]=size(wl);
    g0=0.089;
    g1=0.1245;
    u=(-g0+(g0^2+4*g1.*alphaw_shade).^0.5)/2/g1;
    ki=log10((alphaw_shade(2)+alphaw_shade(3))/(alphaw_shade(4)+5*alphaw_shade(5)*alphaw_shade(5)/alphaw_shade(4)));

    a_ref = aw(id_ref) + 10.^(-1.146 - 1.366 * ki - 0.469*ki.^2);
    bbp550=u(id_ref)*a_ref/(1-u(id_ref));

    ita=2.0*(1-1.2*exp(-0.9*alphaw_shade(2)/alphaw_shade(4)));
    qaa_bbp(1:Wn)=bbp550.*((550./wl(1:Wn)).^ita);
    qaa_a(1:Wn)=(1-u).*(qaa_bbp)./u;
    qaa_a_id=qaa_a(id_ref);

    alphaw_shade0=alphaw_shade(id_ref);
    x0=bbp550;
    LB=0.000001;
    UB=0.5;
    options=optimset('largescale','on','display','iter','tolx',1e-10,'tolfun',1e-12,'Algorithm','active-set','MaxFunEvals',10000,'MaxIter',10000,'Display','off');
    [xx,fval]=fmincon(@(x0) optimization_Ew_pre(x0,qaa_a_id,bbw(id_ref),solz,alphaw_shade0),x0,[],[],[],[],LB,UB,[],options);
    bbp550=xx;

%% End QAA application
%%
%## initial guess implemented###     
%## five variables: Aph440, Adg440, Bbp555, sdg, eta

% bbp550= 0.01; 
P = 0.1;
% G = P; 
G=qaa_a_id;
Sdg = 0.015;
eta=1;

x0=[P,G,bbp550,Sdg,eta];  %initial guess of four variables   
%## spectral optimization ### 
options=optimset('largescale','on','display','iter','tolx',1e-10,'tolfun',1e-12,'Algorithm','active-set','MaxFunEvals',10000,'MaxIter',10000,'Display','off');
% LB=[1e-5;1e-5;1e-5;0.008;-0.5];             %lower boundary
% UB=[5;5;1;0.02;2.5];                         %upper boundary  
LB=[1e-5;1e-5;max(bbp550*0.5,1e-7);0.008;-0.5]; %lower boundary
UB=[5;5;min(bbp550+0.002,bbp550*2.0);0.02;2.5];         %upper boundary  

[xx,~]=fmincon(@(x0) shd_error_opt(x0,wl,aw,bbw,solz,a0,a1,alw_shade),x0,[],[],[],[],LB,UB,[],options); 
[alw_corr,~,shd_error,a,bb,aph,ag]=model_alphaw(xx,wl,aw,bbw,solz,a0,a1);
 

end