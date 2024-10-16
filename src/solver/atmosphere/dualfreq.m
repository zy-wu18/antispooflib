function [rhos_IF, dIons, obs_indice] = dualfreq(obs, uobs, eph_dict)
% Calculate ionospheric delay using dual frequency model
% args  :   3xM double  ps      [m], satellite ECEF position [x, y, z]
%           3XM double  vs      [m/s], satellite ECEF velocity [x, y, z]
%           1xM double  dts     [s], satellite clock fix
%           1xM double  drhos   [], satellite pseodorange rates
%           1xM double  cnrs    [dB·Hz], satellite Carrier-Noise-Ratios
%           1x2M struct obs     satellites dual frequency observation
%           dictionary  eph_dict  ephemeris data struct
% return:   Mx1 double  dIons       [s], ionospheric delay
%           1xM double  rhos_IF     [m], satellite pseodoranges
%           1xM struct  obs_indice  indice of DF satellites observation used

%delete satellites which have no L2C signal
obs_f1 = obs(1:2:end);
obs_f2 = obs(2:2:end);
obs_f1 = obs_f1(~isnan([obs_f2.Rho]));
if (length(uobs) > length(obs_f1))
    obs_indice = ~isnan([obs_f2.Rho]);
else
    obs_indice = 1:length(uobs);
end
obs_f2 = obs_f2(obs_indice);

c = 2.99792458e8;
rhos_f1 = [obs_f1.Rho];
rhos_f2 = [obs_f2.Rho];
M = sum(obs_indice); 
dIons = zeros(2, M);
rhos_IF = zeros(1, M);
%calculate Ion delay
for m = 1:M
    opt = uobs(m).Sys;
    switch opt
        case 'G'
            gamma = (obs_f1(m).Fc/obs_f2(m).Fc)^2;
            rhos_IF(m) = (rhos_f2(m) - gamma * rhos_f1(m)) / (1 - gamma);
            dIons(1, m) = (rhos_f2(m) - rhos_f1(m)) / (gamma - 1);
            dIons(2, m) = gamma* (rhos_f2(m) - rhos_f1(m)) / (gamma - 1);
        case 'C'
            TGD = eph_dict(sprintf("%c%02d", uobs(m).Sys, uobs(m).PRN)).TGD;
            gamma = (obs_f1(m).Fc/obs_f2(m).Fc)^2;
            rhos_f1(m) = rhos_f1(m) - c*TGD(1);
            rhos_IF(m) = (rhos_f2(m) - gamma * rhos_f1(m)) / (1 - gamma);
            dIons(1, m) = (rhos_f2(m) - rhos_f1(m)) / (gamma - 1);
            dIons(2, m) = gamma* (rhos_f2(m) - rhos_f1(m)) / (gamma - 1);
        otherwise
    end
end
dIons = dIons(:, 1) / c;
