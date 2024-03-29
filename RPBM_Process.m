function [RPBM,uRPBM]=RPBM_Process(X,C)
% Returns a structure containing relevant biophysical parameters derived
% from RPBM fit.
%
% (c) Gregory Lemberskiy 2018 gregorylemberskiy@gmail.com


RPBM.D0         = X(1);                                                    %Free Diffusivity [um2/ms]
RPBM.tau        = X(2);                                                    %Tau [ms]from Fit
RPBM.zeta       = X(3);                                                    %Zeta from Fit

RPBM.tortuosity = 1+RPBM.zeta;                                             %Tortuosity
RPBM.L          = sqrt(RPBM.D0*RPBM.tau);                                  %Diffusion Length [um]
RPBM.a          = 2*RPBM.L/RPBM.zeta;                                      %Diameter [um]
RPBM.kappa      = RPBM.D0/(2*RPBM.L);                                      %Permeability [um/ms]
RPBM.SV         = RPBM.zeta*2/RPBM.L;                                      %Surface-to-Volume Ratio [1/um]
RPBM.TD         = RPBM.a^2/(2*RPBM.D0);                                    %Dwell Time [ms]
RPBM.TR         = RPBM.a/(2*RPBM.kappa);                                   %Residence Time [ms]

%% Error Bounds

if C(1,1)==C(1,2)
    uRPBM.D0      = C(1,1);
else
    uRPBM.D0      = (C(1,2)-RPBM.D0 )/2;
end

uRPBM.tau         = (C(2,2)-RPBM.tau)/2;
uRPBM.zeta        = (C(3,2)-RPBM.zeta)/2;

uRPBM.tortuosity = uRPBM.zeta;

uRPBM.L           = 0.5*sqrt((uRPBM.D0/RPBM.D0)^2 ...
        +(uRPBM.tau/RPBM.tau)^2)*(RPBM.D0*RPBM.tau);
uRPBM.a           = 2  *sqrt((uRPBM.L/RPBM.L)^2+(uRPBM.zeta/RPBM.zeta)^2);
uRPBM.kappa       = 0.5*sqrt((uRPBM.L/RPBM.L)^2+(uRPBM.D0/RPBM.D0)^2);
uRPBM.SV          = 0.5*sqrt((uRPBM.L/RPBM.L)^2+(uRPBM.zeta/RPBM.zeta)^2);
uRPBM.TD          = 0.5*sqrt((2*uRPBM.a/RPBM.a)^2+(uRPBM.D0/RPBM.D0)^2);
uRPBM.TR          = 0.5*sqrt((uRPBM.a/RPBM.a)^2+(uRPBM.kappa/RPBM.kappa)^2);

stringRPBM(RPBM,uRPBM)

end


function stringRPBM(RPBM,uRPBM)
%% Prints various RPBM parameters into the command window
%

F = fieldnames(RPBM);
STRP=[];
S=sprintf(['-----  Fitted Results  -----']);

for i=1:length(F)
    STRP=strcat(STRP,sprintf('\n %s = %3.2f +/- %3.2f\n',F{i},RPBM.(F{i}),uRPBM.(F{i})));
end

disp([S,STRP])

end

