clear all; close all; clc;
%% Sample RPBM
% Sample fitting and diffusion parameter estimation from the random 
% permeable barrier model (RPBM) on D(t). This sample script calls
% get_Dt_RPBM to estimate parameters tau and zeta. Afterwards, RPBM_process
% is used to parse these parameters and output relevant biophysical. 
%
%
% Novikov   et al. Nature Physics 7:508 (2011)
% Fiermeans et al. NMR in Biomed 30:e3612 (2017)
%
% (c) Gregory Lemberskiy 2017 gregorylemberskiy@gmail.com


%% Input - Sample Soleus Muscle D(t) following DTI Tensor Fitting
% Soleus ROI from a DWI acquisition using 3x3x5 mm^3 voxels, 20 directions
% at b=500 s/mm^2 and 3 b=0. Outlier rejection was performed using IRESTORE
% and tensor estimation was performed using weighted linear least squares.
%
% Chang   et al. Magn Reson Med 68: 1654?1663 (2012) 
% Veraart et al. Neuroimage     81:  335-46   (2013)

time = [57,75,100,200,350,500,750,1000,1200];                              %Diffusion Time
DL   = [1.7758,1.7635,1.7511,1.7189,1.7110,1.6791,1.6780,1.6580,1.6864];   %Longitudinal Diffusivity
DR   = [1.3523,1.2932,1.2635,1.1738,1.1142,1.0645,1.0133,0.9818,1.0002];   %Radial Diffusivity


%% Fitting by Fixing D0 to L1(t~inf) 
% In principle, D0 should be equal to L1 at long times. However, in the 
% presense of noisy data this the eigenvalues will repel from one another
% causing a systematic bias on parameter estimation. 
%
% Conversely if the structure of interest is highly restricted (small 
% diameters), diffusion times may be too long to accurately estimate D0. 

Dfix = mean(DL(time>200));

fun_fix = @(F,xdata) Dfix*get_Dt_RPBM(xdata./F(1),F(2));
lb = [0  0]; ub = [Inf Inf]; x0 = [100 2];

[Xfix,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(fun_fix,x0,time,DR,lb,ub);
ci_fix = nlparci(Xfix,resid,'jacobian',J);

[RPBM_fix,uRPBM_fix]=RPBM_Process([Dfix,Xfix],[0 0;ci_fix]);


%% Fitting by Varying D0
% If there is sufficient SNR and sampling of t, fitting over 3 parameters 
% would be the ideal approach. Given your acquisition, I recommend trying
% both approaches [fitting/fixing D0] to see if they are in agreement. 

fun_vary = @(F,xdata) F(1)*get_Dt_RPBM(xdata./F(2),F(3));
lb  =  [0      0   0]; ub  =  [3    Inf Inf]; x0  =  [Dfix 100   2];

[Xvary,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(fun_vary,x0,time,DR,lb,ub);
ci_vary = nlparci(Xvary,resid,'jacobian',J);

[RPBM_vary,uRPBM_vary]=RPBM_Process(Xvary,ci_vary);


%% Plotting Results

txn=linspace(0,14000,20000);
VFIX=fun_fix(Xfix,txn);
VFIT=fun_vary(Xvary,txn);
%%
figure; 
subplot(1,2,1)
plot(time,DL,'rs','markerfacecolor','r','markersize',10); hold on;
plot(time,DR,'ro','markerfacecolor','r','markersize',10);
plot(txn,VFIX,'r--','linewidth',2)
plot(txn,VFIT,'b--','linewidth',2)
xlabel('Time [ms]','fontsize',20)
ylabel('D_{||}, D_\perp','fontsize',20)
set(gca,'linewidth',2,'fontweight','bold')
pbaspect([1 1 1])
ylim([0.9 1.9])
xlim([0 1400])

subplot(1,2,2)
plot(time,DL,'rs','markerfacecolor','r','markersize',10); hold on;
plot(time,DR,'ro','markerfacecolor','r','markersize',10);
plot(txn,VFIX,'r--','linewidth',2)
plot(txn,VFIT,'b--','linewidth',2)
xlabel('Time [ms]','fontsize',20)
ylabel('D_{||}, D_\perp','fontsize',20)
set(gca,'linewidth',2,'fontweight','bold')
set(gca,'xscale','log')
set(gca,'xtick',[1 10 100 500 2500 10000])
title('Semilog X','fontsize',20)
pbaspect([1 1 1])
ylim([0.9 1.9])
xlim([txn(1) txn(end)])
