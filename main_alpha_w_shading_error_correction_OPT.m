clear
clc;
load('MC007.mat');
% load('MC016.mat');
% load('MC037.mat');
wl=MC.wl;
alw_shd=MC.alw_shade;
sza=30; % MC.solz;
sno=size(alw_shd,1);

% load("G_LUT.mat");
%### start optimization for each spectra ####
for i=1:sno
    disp(['Current processing ' num2str(i) '/' num2str(sno) ' Sample!']);
    solz=sza;     
    alw_shade=alw_shd(i,:); 
    alw_shade(alw_shade<1e-6)=1e-6;
%     [alw_corr_temp,shd_error_temp,a_temp,bb_temp,~,~]=shade_corr_alw(wl,alw_shade,solz,G_LUT);
    [alw_corr_temp,shd_error_temp,a_temp,bb_temp,~,~]=shade_corr_alw(wl,alw_shade,solz);
    alw_corr(i,:)=alw_corr_temp;    
    shd_error(i,:)=shd_error_temp;
    a_der(i,:)=a_temp;
    bb_der(i,:)=bb_temp;
    
end

% 
MC.alw_corr=alw_corr;
MC.shd_error=shd_error;
MC.a_der=a_der;
MC.bb_der=bb_der;

s1=mfilename('fullpath')
[filepath,~,~]=fileparts(s1)
cd (filepath)
save('MC007_corr_ZH.mat','MC');
% save('MC037_corr_ZH.mat','MC');
% save('MC016_corr_ZH.mat','MC');
% save('MC037_corr.mat','MC');