function [RPBM,uRPBM]=RPBM_Process(X,C)
% Returns a structure containing relevant biophysical parameters derived
% from RPBM fit.
%
% (c) Gregory Lemberskiy 2018 gregorylemberskiy@gmail.com

RPBM.D0         = X(1);                          %Free Diffusivity [um2/ms]
RPBM.tau        = X(2);                                   %Tau [ms]from Fit
RPBM.zeta       = X(3);                                      %Zeta from Fit

RPBM.tortuosity = 1+RPBM.zeta;                                  %Tortuosity
RPBM.L          = sqrt(RPBM.D0*RPBM.tau);            %Diffusion Length [um]
RPBM.a          = 2*RPBM.L/RPBM.zeta;                      %Diameter [um]
RPBM.kappa      = RPBM.D0/(2*RPBM.L);               %Permeability [um/ms]
RPBM.SV         = RPBM.zeta*2/RPBM.L;     %Surface-to-Volume Ratio [1/um]
RPBM.TD         = RPBM.a^2/(2*RPBM.D0);                    %DwL Time [ms]
RPBM.TR         = RPBM.a/(2*RPBM.kappa);

%% Error Bounds

if C(1,2)==0
    uRPBM.D0         = 0;
else
    uRPBM.D0         = (C(1,2)-RPBM.D0)/2;
end

uRPBM.tau        = (C(2,2)-RPBM.tau)/2;
uRPBM.zeta       = (C(3,2)-RPBM.zeta)/2;

uRPBM.tortuosity = uRPBM.zeta;
if C(1,2)==0
    uRPBM.L          = 0.5*RPBM.D0*uRPBM.tau/RPBM.tau;
else
    uRPBM.L          = 0.5*sqrt((uRPBM.D0/RPBM.D0)^2+(uRPBM.tau/RPBM.tau)^2)*(RPBM.D0*RPBM.tau);
end 
uRPBM.a          = 2  *sqrt((uRPBM.L/RPBM.L)^2+(uRPBM.zeta/RPBM.zeta)^2);
uRPBM.kappa      = 0.5*sqrt((uRPBM.L/RPBM.L)^2+(uRPBM.D0/RPBM.D0)^2);
uRPBM.SV         = 0.5*sqrt((uRPBM.L/RPBM.L)^2+(uRPBM.zeta/RPBM.zeta)^2);
uRPBM.TD         = 0.5*sqrt((2*uRPBM.a/RPBM.a)^2+(uRPBM.D0/RPBM.D0)^2);
uRPBM.TR         = 0.5*sqrt((uRPBM.a/RPBM.a)^2+(uRPBM.kappa/RPBM.kappa)^2);


end