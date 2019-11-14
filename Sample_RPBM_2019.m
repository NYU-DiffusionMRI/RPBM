clear all; close all; clc;
%% Sample RPBM
% Sample fitting and diffusion parameter estimation from the random 
% permeable barrier model (RPBM) on D(t). This sample script calls
% get_Dt_RPBM to estimate parameters tau and zeta. Afterwards, RPBM_process
% is used to parse these parameters and output relevant biophysical. 
%
% Novikov   et al. Nature Physics 7:508 (2011)
% Fiermeans et al. NMR in Biomed 30:e3612 (2017)
%
% (c) Gregory Lemberskiy 2019 gregorylemberskiy@gmail.com

%% Input - Sample Soleus Muscle D(t) following DTI Tensor Fitting
% Soleus ROI from a DWI acquisition using 3x3x5 mm^3 voxels, 20 directions
% at b=500 s/mm^2 and 3 b=0. Denoising, degibbs, and outlier rejection was 
% performed using IRESTORE and tensor estimation was performed using 
% weighted linear least squares.
% 
% Chang   et al. IRESTORE. Magn Reson Med   68: 1654-1663 (2012) 
% Veraart et al. WLLS    . Neuroimage       81:  335-46   (2013)
% Veraart et al. Denoising.Neuroimage      142:  394-406  (2016)
% Kellner et al. DeGibbs.  Magn Reson Med   76: 1574?1581 (2016)

time = [57,75,100,200,350,500,750,1000,1200];                              %Diffusion Time
DL   = [1.7758,1.7635,1.7511,1.7189,1.7110,1.6791,1.6780,1.6580,1.6864];   %Longitudinal Diffusivity
DR   = [1.3523,1.2932,1.2635,1.1738,1.1142,1.0645,1.0133,0.9818,1.0002];   %Radial Diffusivity

uDL=[0.1438,0.1579,0.1407,0.1542,0.1686,0.1596,0.1886,0.2093,0.2460];
uDR=[0.1050,0.0902,0.0858,0.0993,0.1075,0.1133,0.1348,0.1862,0.2045];

I.t                =time;
I.L1               =DL;
I.DR               =DR;
I.dictionarysize   =1000;
I.dictionarynumber =500;
I.SNR              =35;

Dfix =mean(DL(time>200),2);
sDfix= std(DL(time>200),[],2);

%% Dictionary Generation
%[sig,Z,T]=GS_RPBM_dictionary(I);
load('/Users/gl744/Documents/Scripts/Processing/RPBM_DICTIONARY.mat');
load('Img.mat')

SNR=35;
Nx=size(DR,1);
for nx=1:Nx
sig(:,nx,:,:)=Dfix(nx).*sig;
end

%Adding noise approximately
sigma   = 3/SNR;
noise_1 = normrnd(0, sigma, size(sig));
sig_n   = sig + noise_1;

%% Dictionary Matching
clear ZZ TT
dsize=size(sig,3);
dnum =size(sig,4);

for ic=1:dnum
    for nx=1:Nx
        [~,ix]  = min(sum(abs(sig_n(:,:,:,ic)-DR(nx,:)').^2,1));
        %Outlier Detection based on standard error
        A=abs((sig_n(:,:,ix,ic)-DR(nx,:)')./DR(nx,:)');
        O=(single(isnan(A))+single(A>25))>1.0;
        if sum(O)>0
            [~,ix]  = min(sum((abs(sig_n(~O,:,:,ic)-DR(nx,~O)')).^2,1));
        end
        ZZ(nx,ic)  = Z( ix,ic);
        TT(nx,ic)  = T( ix,ic);
    end
end

X=[Dfix,median(TT),median(ZZ)];
C=[Dfix-sDfix Dfix+sDfix;median(TT)-std(TT) median(TT)+std(TT);median(ZZ)-std(ZZ) median(ZZ)+std(ZZ);];
[RPBM,uRPBM]=RPBM_Process(X,C);

%% Plotting D(t)
txn=exp(linspace(log(1),log(100000),1000));
RPBM_FIT = RPBM.D0.*get_Dt_2015(txn./ RPBM.tau, RPBM.zeta);
uRPBM_FIT=uRPBM.D0.*get_Dt_2015(txn./uRPBM.tau,uRPBM.zeta);

figure(1)
subplot(1,1,1)
plot(time,DL,'rs','markerfacecolor','r','markersize',10); hold on;
plot(time,DR,'ro','markerfacecolor','r','markersize',10);
x1(1)=plot(txn,RPBM_FIT,'r--','linewidth',2);
plot(txn,RPBM_FIT+uRPBM_FIT,'-','color',[1,0,0].*0.5);                     %Error Lines
plot(txn,RPBM_FIT-uRPBM_FIT,'-','color',[1,0,0].*0.5);                     %Error Lines
set(gca,'linewidth',1,'ticklabelinterpreter','latex','fontsize',15)
xlabel('$t\,[ms]$','fontsize',25,'interpreter','latex')
ylabel('$D_{||},\,D_\perp$','fontsize',25,'interpreter','latex')
pbaspect([1 1 1])
ylim([0.925    1.825]);
xlim([0     1400    ]);

%% Dictionary Matching Full Image 
% [~5 minutes on Macbook Pro 15" 2013], 
%  Feel free to skip this step and load the result in the next section. 
clear ZZ_ TT_
dsize=size(sig,3);
dnum =size(sig,4);
[Nx,Ny,Nz,Nt]=size(DRimg);
L1img_=reshape(L1img,[Nx*Ny*Nz,Nt]);
DRimg_=reshape(DRimg,[Nx*Ny*Nz,Nt]);

D0=mean(L1img_(:,time>200),2);
uD0=std(L1img_(:,time>200),[],2);
Nm=Nx*Ny*Nz;
SNR=35;

tic
parfor nx=1:Nm
    nx
    sigma   = 3/SNR;
    noise_1 = normrnd(0, sigma, size(sig));
    sig_n=D0(nx).*sig+ noise_1;
    %Adding noise approximately
    for ic=1:dnum
        if nnz(DRimg_(nx,:))>0
            [~,ix]  = min(sum(abs(sig_n(:,:,:,ic)-DRimg_(nx,:)').^2,1));
            %Outlier Detection based on standard error
            A=abs((sig_n(:,:,ix,ic)-DRimg_(nx,:)')./DRimg_(nx,:)');
            O=(single(isnan(A))+single(A>25))>1.0;
            if sum(O)>0
                [~,ix]  = min(sum((abs(sig_n(~O,:,:,ic)-DRimg_(nx,~O)')).^2,1));
            end
            %DD0(ib) = D0_(ix);
            ZZ_(nx,ic) = Z( ix,ic);
            TT_(nx,ic) = T( ix,ic);
        else
            ZZ_(nx,ic) = 0;
            TT_(nx,ic) = 0;
        end
    end
end
toc
%% Plotting RPBM FIT
load('RPBMMatch.mat')
clear RPBMimg uRPBMimg a_img kappa_img SV_img zeta_img tau_img TD_img TR_img tortuosity_img
mTT=median(TT_   ,2);
sTT=   std(TT_,[],2);
mZZ=median(ZZ_   ,2);
sZZ=   std(ZZ_,[],2);

for nm=1:Nm
    X=[D0(nm),mTT(nm),mZZ(nm)];
    C=[D0(nm)-uD0(nm) D0(nm)+uD0(nm);mTT(nm)-sTT(nm) mTT(nm)+sTT(nm);mZZ(nm)-sZZ(nm) mZZ(nm)+sZZ(nm);];
    [RPBMimg{nm},uRPBMimg{nm}]=RPBM_Process(X,C);
    a_img(nm)            =RPBMimg{nm}.a;
    kappa_img(nm)        =RPBMimg{nm}.kappa;
    zeta_img(nm)         =RPBMimg{nm}.zeta;
    tau_img(nm)          =RPBMimg{nm}.tau;
    SV_img(nm)           =RPBMimg{nm}.SV;
    kappa_img(nm)        =RPBMimg{nm}.kappa;
    TD_img(nm)           =RPBMimg{nm}.TD;
    TR_img(nm)           =RPBMimg{nm}.TR;
    tortuosity_img(nm)   =RPBMimg{nm}.tortuosity;
end

PM(:,:,1)=reshape(a_img,[Nx,Ny,Nz]);
PM(:,:,2)=reshape(kappa_img,[Nx,Ny,Nz]);
PM(:,:,3)=reshape(D0,[Nx,Ny,Nz]);
PM(:,:,4)=reshape(zeta_img,[Nx,Ny,Nz]);
PM(:,:,5)=reshape(tau_img,[Nx,Ny,Nz]);
PM(:,:,6)=reshape(TD_img,[Nx,Ny,Nz]);
PM(:,:,7)=reshape(TR_img,[Nx,Ny,Nz]);
PM(:,:,8)=reshape(tortuosity_img,[Nx,Ny,Nz]);

XL=[0 40;
    0 0.15;
    0 3;
    0 6;
    0 1000;
    0 200
    0 600
    0 10];
XLname={'$\bar{a}\,[\mu m]$','$\kappa\,[\mu m/ms]$','$D_0\,[\mu m^2/ms]$','$\zeta$','$\tau\,[ms]$','$T_{D}\,[ms]$','$T_{R}\,[ms]$','Tortuosity'};
figure(5);
for pm=1:8
    subplot(2,4,pm)
    imagesc(rot90(PM(:,:,pm)),XL(pm,:))
    axis('equal');
    xlim([1 size(PM,1)]);
    ylim([1 size(PM,2)]);
    set(gca,'xtick',[],'ytick',[]);
    title([XLname(pm)],'interpreter','latex','fontsize',25);
    colormap hot
    colorbar('location','southoutside','fontsize',10,'linewidth',1)
end
set(gcf,'Position',[50 900 1000 500])
