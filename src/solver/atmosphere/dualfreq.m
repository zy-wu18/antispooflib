function [rhosIF, dions, uobs] = dualfreq(obs, uobs, eph_dict)
% Calculate ionospheric delay using dual frequency model
% args  :   1x2M obs_t  obs         observation, both 2nd frequency and
%                                   filtered ones included
%           1xM  obs_t  uobs        used observation, only main(1st) ones 
%           dictionary  eph_dict    ephemeris data struct
% return:   Mx1 double  dIons       [s], ionospheric delay
%           1xM double  rhos_IF     [m], satellite pseodoranges

obs1 = uobs;
if(isempty(obs1))
    obs2 = obs1;
else
    obs2 = [];
end

p = 1;
for m = 1:length(uobs)
    key_m = sprintf('%c%02d', obs1(m).Sys, obs1(m).PRN);
    while(p < length(obs))
        key_p = sprintf('%c%02d', obs(p).Sys, obs(p).PRN);
        if(strcmp(key_m, key_p))
            p = p + 1;
            obs2 = [obs2, obs(p)]; %#ok<AGROW> 
            assert(obs(p).Fc ~= obs1(m).Fc);
            break;
        end
        p = p + 1;
    end
end
assert(length(obs1) == length(obs2));

c = 2.99792458e8;
rhos1 = [obs1.Rho];
rhos2 = [obs2.Rho];
M = length(obs1);
dions = zeros(2, M);
rhosIF = zeros(1, M);
%calculate Ion delay
for m = 1:M
    sys = uobs(m).Sys;
    switch sys
        case 'G'
            gamma = (obs1(m).Fc/obs2(m).Fc)^2;
            rhosIF(m) = (rhos2(m) - gamma * rhos1(m)) / (1 - gamma);
            dions(1, m) = (rhos2(m) - rhos1(m)) / (gamma - 1);
            dions(2, m) = gamma* (rhos2(m) - rhos1(m)) / (gamma - 1);
        case 'C'
            TGD = eph_dict(sprintf("%c%02d", uobs(m).Sys, uobs(m).PRN)).TGD;
            gamma = (obs1(m).Fc/obs2(m).Fc)^2;
            rhos1(m) = rhos1(m) - c*TGD(1);
            rhosIF(m) = (rhos2(m) - gamma * rhos1(m)) / (1 - gamma);
            dions(1, m) = (rhos2(m) - rhos1(m)) / (gamma - 1);
            dions(2, m) = gamma* (rhos2(m) - rhos1(m)) / (gamma - 1);
        otherwise
            % TODO
    end
end
dions = dions(1,:) / c; % default: fix on L1 pseudo range
