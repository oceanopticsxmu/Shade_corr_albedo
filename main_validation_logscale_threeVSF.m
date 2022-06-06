clear;
clc;

ti=['a';'b';'c';'d';'e';'f'];
close all;
figure;
set(gcf,'position',[100,400, 950,600]);
mc=[0.85,0.33,0.10;0.00,0.45,0.74;0.93,0.69,0.13;0,0.145,1];


vsf='007';
load(['MC' vsf '_corr_ZH.mat']);
alw1=MC.alw;
alw_cor1=MC.alw_corr;
wl=MC.wl;
sno1=size(alw1,1);

vsf='016';
load(['MC' vsf '_corr_ZH.mat']);
alw2=MC.alw;
alw_cor2=MC.alw_corr;   
sno2=size(alw2,1);

vsf='037';
load(['MC' vsf '_corr_ZH.mat']);
% load(['MC' vsf '_corr_backup.mat']);
alw3=MC.alw;
alw_cor3=MC.alw_corr;
sno3=size(alw3,1);
alw=[alw1;alw2;alw3];
alw_cor=[alw_cor1;alw_cor2;alw_cor3];

vsfidx=[1+zeros(sno1,1),2+zeros(sno2,1),3+zeros(sno3,1)];
sno=sno1+sno2+sno3;
idx=1:1:sno;

for i= 1:length(wl)  
    idx_temp=find(alw(:,i)>0);
    idx=intersect(idx,idx_temp);     
end
% load('QAAvalididx.mat');

for sn=1:length(wl)
    
    subplot(2,3,sn);  
    yy=alw_cor(:,sn);    
%     idx2=find(yy>0 & yy<1);
%     idx=intersect(idx,idx2);
 
    x=alw(idx,sn);
    y=yy(idx);         
    vidx=vsfidx(idx);
    vidx1=find(vidx==1);
    vidx2=find(vidx==2);
    vidx3=find(vidx==3);
    
    plot(x(vidx1),y(vidx1),'o','color',mc(1,:));
    hold on;
    plot(x(vidx2),y(vidx2),'s','color',mc(2,:));
    plot(x(vidx3),y(vidx3),'+','color',mc(3,:));
    box on;
    set(gca,'xscale','log','yscale','log','fontname','times new roman','fontsize',13);  
    
    uplim=ceil(log10(max(max([x,y]))));
    lowlim=floor(log10(min(min([x,y]))));
    xlim([10^lowlim 10^uplim]);
    ylim([10^lowlim 10^uplim]);
    h=refline(1,0); set(h,'linestyle','--','color',[0.1 0.1 0.1],'linewidth',1.5);
    grid on;
%   set(gca,'FontSize',13,'FontName','Times New Roman');
    
    xlabel('Known \alpha_w','FontSize',12,'FontName','Times New Roman','fontweight','bold');
    ylabel('Corrected \alpha_w','FontSize',12,'FontName','Times New Roman','fontweight','bold');
    title(['(' ti(sn) ') \lambda = ' num2str(wl(sn)) ' nm'],'FontSize',12,'FontName','Times New Roman','fontweight','bold');
    
    xmax=uplim;
    xmin=lowlim;     
    
    step=abs(uplim-lowlim)./100;
    xshift=95;
    yshift=21; 
    [m,b,rd,sm,sb]=lsqfitgm(x,y); 
    R=ceil((rd.^2)*1000)./1000;
    MPD=100*median(abs(y./x-1));
%   RMSE=sqrt(mean((y-x).^2));    
    bias=100*median((y./x-1));
        num = round(m,2);
        str=['Slope = ' num2str(num)];
        text(10^(xmax-xshift*step),10^(xmin+4.4*yshift*step),str,'FontSize',12,'FontName','Times New Roman');
   
        num = round(R,2);
        str=['{\itR}^2 = ' num2str(num)];
        text(10^(xmax-xshift*step),10^(xmin+3.9*yshift*step),str,'FontSize',12,'FontName','Times New Roman');
   
        num = round(MPD,1);
        str=['MAPD = ' num2str(num) '%'];
        text(10^(xmax-xshift*step),10^(xmin+3.4*yshift*step),str,'FontSize',12,'FontName','Times New Roman'); 
        
%         num = round(bias,1);
%         str=['{\it{bias}} = ' num2str(num) '%'];
%         text(10^(xmax-xshift*step),10^(xmin+2.9*yshift*step),str,'FontSize',12,'FontName','Times New Roman'); 
%         
        num = round(length(idx),1);
        str=['{\itN} = ' num2str(num)];
        text(10^(xmax-xshift*step),10^(xmin+2.9*yshift*step),str,'FontSize',12,'FontName','Times New Roman'); 
        
        if sn==1             
            lgtxt={'{\itb_{bpr}007}';'{\itb_{bpr}016}';'{\itb_{bpr}037}'};
            legend(lgtxt,'FontSize',12,'FontName','Times New Roman','box','off','Location','southeast');
        end
        
end

% print(gcf,'-r250','-dbitmap','Val_shd_cor_3vsf_QAAidx.bmp');
% print(gcf,'-r250','-dbitmap','Val_shd_cor_3vsf_QAAidx_backup.bmp');
% close all;


% print(gcf,'-r250','-dbitmap',['Val_alw_corr_QAA_' vsf '.bmp']);