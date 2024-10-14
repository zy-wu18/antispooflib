function [rhos_IF, drhos_IF, ps_IF, vs_IF, dts_IF, cnrs_IF, dIons, uobs_out] = ...
    DoubleFreq(ps, vs, dts, drhos, cnrs, obs, uobs, opt)
%delete satellites which have no L2C signal
obs_f1 = obs(1:2:end);
obs_f2 = obs(2:2:end);

obs_f1 = obs_f1(~isnan([obs_f2.Rho]));
ps_IF = ps(:, ~isnan([obs_f2.Rho]));
vs_IF = vs(:, ~isnan([obs_f2.Rho]));
dts_IF = dts(~isnan([obs_f2.Rho]));
drhos_IF = drhos(~isnan([obs_f2.Rho]));
if (length(uobs) > length(obs_f1))
    uobs_out = uobs(~isnan([obs_f2.Rho]));
else
    uobs_out = uobs;
end
cnrs_IF = cnrs(~isnan([obs_f2.Rho]));
obs_f2 = obs_f2(~isnan([obs_f2.Rho]));

c = 2.99792458e8; % [m/s], speed of light
rhos_f1 = [obs_f1.Rho] ;
rhos_f2 = [obs_f2.Rho] ;

prn_num = length(obs_f1); 
dIons = zeros(32, 2);
rhos_IF = zeros(prn_num,1);

%calculate Ion delay
switch opt
    case 'IonoFree'
        for k = 1:prn_num 
            gamma = (obs_f1(k).Fc/obs_f2(k).Fc)^2;
            rhos_IF(k) = (rhos_f2(k) - gamma * rhos_f1(k)) / (1 - gamma);
            dIons(k, 1) = (rhos_f2(k) - rhos_f1(k)) / (gamma - 1);
            dIons(k, 2) = gamma* (rhos_f2(k) - rhos_f1(k)) / (gamma - 1);
        end
    case 'IonoFree_CarrPhase'
        c = 2.99792458e8;
        for k = 1:prn_num 
            lambda_f1 = c/obs_f1(k).Fc;
            lambda_f2 = c/obs_f2(k).Fc;

            gamma = (obs_f1(k).Fc)^2/((obs_f1(k).Fc)^2 - (obs_f2(k).Fc)^2);
            rhos_IF(k) = gamma * lambda_f1*obs_f1(k).AcPh + (1 - gamma)*lambda_f2*obs_f2(k).AcPh;
            dIons(k, 1) = (1 - gamma)*(lambda_f1*obs_f1(k).AcPh - lambda_f2*obs_f2(k).AcPh);
            dIons(k, 2) = gamma*(lambda_f2*obs_f2(k).AcPh -lambda_f1*obs_f1(k).AcPh);
           
        end
    otherwise
end
dIons = dIons(:, 1) / c;
rhos_IF = rhos_IF';