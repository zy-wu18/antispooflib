function [rhos_IF, drhos_IF, ps_IF, vs_IF, dts_IF, cnrs_IF, dIons, uobs_out] = ...
    dualfreq(ps, vs, dts, drhos, cnrs, obs, uobs, eph_dict)
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
%           1xM double  drhos_IF    [], satellite pseodorange rates
%           3xM double  ps_IF       [m], satellite ECEF position [x, y, z]
%           3XM double  vs_IF       [m/s], satellite ECEF velocity [x, y, z]
%           1xM double  dts_IF      [s], satellite clock fix%           
%           1xM double  cnrs_IF     [dB·Hz], satellite Carrier-Noise-Ratios
%           1xM struct  uobs_out    satellites observation used in PVT calculation

%delete satellites which have no L2C signal
obs_f1 = obs(1:2:end);
obs_f2 = obs(2:2:end);
obs_f1 = obs_f1(~isnan([obs_f2.Rho]));
if (length(uobs) > length(obs_f1))
    uobs_out = uobs(~isnan([obs_f2.Rho]));
    ps_IF = ps(:, ~isnan([obs_f2.Rho]));
    vs_IF = vs(:, ~isnan([obs_f2.Rho]));
    dts_IF = dts(~isnan([obs_f2.Rho]));
    drhos_IF = drhos(~isnan([obs_f2.Rho]));
    cnrs_IF = cnrs(~isnan([obs_f2.Rho]));
else
    ps_IF = ps;
    vs_IF = vs;
    dts_IF = dts;
    drhos_IF = drhos;
    cnrs_IF = cnrs;
    uobs_out = uobs;
end
obs_f2 = obs_f2(~isnan([obs_f2.Rho]));

c = 2.99792458e8;
rhos_f1 = [obs_f1.Rho] ;
rhos_f2 = [obs_f2.Rho] ;
prn_num = length(uobs_out); 
dIons = zeros(32, 2);
rhos_IF = zeros(prn_num,1);
%calculate Ion delay
for k = 1:prn_num
    opt = uobs(k).Sys;
    switch opt
        case 'G'
            gamma = (obs_f1(k).Fc/obs_f2(k).Fc)^2;
            rhos_IF(k) = (rhos_f2(k) - gamma * rhos_f1(k)) / (1 - gamma);
            dIons(k, 1) = (rhos_f2(k) - rhos_f1(k)) / (gamma - 1);
            dIons(k, 2) = gamma* (rhos_f2(k) - rhos_f1(k)) / (gamma - 1);
        case 'C'
            TGD = eph_dict(sprintf("%c%02d", uobs(k).Sys, uobs(k).PRN)).TGD;
            gamma = (obs_f1(k).Fc/obs_f2(k).Fc)^2;
            rhos_f1(k) = rhos_f1(k) - c*TGD(1);
            rhos_IF(k) = (rhos_f2(k) - gamma * rhos_f1(k)) / (1 - gamma);
            dIons(k, 1) = (rhos_f2(k) - rhos_f1(k)) / (gamma - 1);
            dIons(k, 2) = gamma* (rhos_f2(k) - rhos_f1(k)) / (gamma - 1);
        otherwise
    end
end
dIons = dIons(:, 1) / c;
rhos_IF = rhos_IF';