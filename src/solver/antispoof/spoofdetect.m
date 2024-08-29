function sdres = spoofdetect(rhos, drhos, ps, vs, dts, cnrs, keys, mthds)
% Spoof detector wrapper, with all available observables as input
% args  :   1xM double  rho     [m], bias-fixed pseudorange observables
%           1xM double  drho    [m/s], bias-fixed pseudorange-rate observables
%           3xM double  ps      [m], satellite position in ECEF
%           3xM double  vs      [m/s], satellite velocity in ECEF
%           1xM double  cnrs    [dBHz], satellite carrier-to-noise ratios
%           raimcfg_t   cfg     (optional)customized configuration
% return:   sdres_t     sdres   name, alarm, test statistics and threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nmthd = length(mthds);
    sdres_t = struct('name', [], 'alarm', [], 'T', [], 'gamma', []);
    sdres = repmat(sdres_t, [1, nmthd]);
    
    for i = 1:nmthd
        sdmthd = mthds(i);
        sdres(i) = sdmthd.detector(rhos, drhos, ps, vs, dts, cnrs, keys, sdmthd.cfg);
    end
end
