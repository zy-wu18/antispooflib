function [ps, vs, dts, rhos, drhos, cnrs, obs_pvt] = ephposfix(obs_raw, eph_dict, pntcfg, pu_ref)
%EPHPOSFIX ephemeris positioning and measurement fix
% args  :   1xM0 obs_t  obs_raw     a single raw observable frame
%           dictionary  eph_dict    string->struct(neph_t, or geph_t) with
%                                   key=sprintf("%c%02d", sys, prn)
%           pntcfg_t    pntcfg      PNT configuration, include 
%                                   elvMask([deg], default=0), <!TODO>
%                                   constellation(char, default='GRJEC'), 
%                                   cnrMask([dBHz], default=0).
% return:   Mx3 double  ps          [m], fixed satellite ECEF position
%           Mx3 double  vs          [m/s], fixed satellite ECEF velocity
%           Mx1 double  rhos        [m], pseudo-range with satellite clock 
%                                   bias fixed.
%           Mx1 double  drhos       [m/s], pseudo-range rate, which is
%                                   calculated as drho = -fd/fc*c
%           Mx1 double  cnrs        [dBHz], carrier-noise-ratio
%           1xM obs_t   obs_pvt     a single filtered observable frame
%                                   with Az and El[rad] attached
% notes:    A) filter the input raw observables(#=M0) according to pntcfg
%           B) calculate the ps, vs, dts of each filtered satellite(#=M)
%           C) fix the original measurements according to dts
%                   rho <= rho + dts*c, drho <= -1.0*fd/fc*c.

    %% filter observations, M0 observables -> M observables/satellites
    obs_pvt = obs_raw; % observables that passed rough check
    M0 = length(obs_raw);
    M = 0;
    elv_mask_off = isempty(pu_ref) || any(isnan(pu_ref));

    for j = 1:M0
        o = obs_raw(j);
        key = sprintf("%c%02d", o.Sys, o.PRN);
        c1 = isKey(eph_dict, key) && (eph_dict(key).Health==0);
        c2 = ~(isempty(o.Fd) || isempty(o.Rho));
        c3 = ~((j>1) && (o.Sys==obs_raw(j-1).Sys) && (o.PRN==obs_raw(j-1).PRN));
        c4 = any(o.Sys == pntcfg.constellation) && o.CNR > pntcfg.cnrMask;
        c5 = isKey(eph_dict, key) && (elv_mask_off || ...
            (satelaz(pu_ref,eph2pvt(o.ObsTime, eph_dict(key)))/pi*180 >= pntcfg.elvMask));

        if(c1 && c2 && c3 && c4 && c5)
            M = M + 1;
            obs_pvt(M) = o;
        end
    end
    obs_pvt = obs_pvt(1:M);

    %% use ephemeris to calculate the PVT of selected satellites, #satu=M
    ps = zeros(3, M);
    vs = zeros(3, M);
    dts = zeros(1, M);
    rhos = zeros(1, M);
    drhos= zeros(1, M);
    cnrs = zeros(1, M);
    c = 2.99792458e8; % [m/s], speed of light
    
    for j = 1:M
        o = obs_pvt(j);
        key = sprintf("%c%02d", o.Sys, o.PRN);
        eph = eph_dict(key);
        tsv = o.ObsTime - o.Rho/c; % transmit time = tlatch - rho/c
        tsv = tsv - eph2clk(tsv, eph);
        [ps(:, j), vs(:, j), dts(j)] = eph2pvt(tsv, eph);
        rhos(j) = o.Rho;
        drhos(j)= -1.0*o.Fd/o.Fc*c;
        cnrs(j) = o.CNR;
        [obs_pvt(j).El, obs_pvt(j).Az] = satelaz(pu_ref, ps(:,j));
    end
end
