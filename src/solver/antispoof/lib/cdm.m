function sdres = cdm(drho, ps, vs, cfg)
% clock drift monitoring
% args  :   1xM double  drho    [m/s], bias-fixed pseudorange-rate
%           3xM double  ps      [m], satellite position in ECEF
%           3xM double  vs      [m/s], satellite velocity in ECEF
%           draimcfg_t  cfg     (optional)customized configuration
% return:   sdres_t     sdres   name, alarm, test statistics and threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(nargin < 3 || isempty(cfg))
        cfg = struct();
        cfg.ObsRate=10;
        warning("cdm: no cfg input, obs rate set to default(%dHz).", cfg.ObsRate);
        cfg.TAccum=1.0;
        cfg.NumPoll=5;
        cfg.Af2Thresh=1e-9; % clock-drift([s/s]) per observable sample
    end
    sdres = struct("name", 'ClkDrfMontoring', 'alarm', NaN, 'T', NaN, 'gamma', NaN);
    N = round(cfg.TAccum*cfg.ObsRate);   
    assert(isscalar(N) && N >= 3);
    K = cfg.NumPoll;    
    assert(isscalar(K) && K >= 1);

    %% Observable accumulation
    [~, ~, ~, ddtu, ~, ~, ~] = lse7pnt(drho, ps, vs);
    persistent i; % thread counter
    persistent ddtu_arr;
    persistent M_arr;
    if(isempty(i))
        i = 0;
        ddtu_arr = zeros(N, 1);
        M_arr = zeros(1, N);
    end
    ddtu_arr(i+1) = ddtu;
    i = mod(i+1, N);
    M_arr(i+1) = length(drho);
    
    if(mod(i, N)~= 0 || min(M_arr) <= 7)
        return;
    end
    
    %% CDM Core
    Gddt = [ones(1, N); (0:N-1)/cfg.ObsRate]';
    af = (Gddt'*Gddt)\Gddt'*ddtu_arr;
    sdres.T = abs(af(2)); % 1st clock drift deviation, [s/s^2]
    sdres.gamma = cfg.Af2Thresh;
    
    %% Output filtering
    persistent alarm_p;
    if(K ~= length(alarm_p))
        alarm_p = zeros(1, K);
    end
    alarm_t = (sdres.T > sdres.gamma);
    alarm_p = [alarm_t, alarm_p(1:end-1)];
    sdres.alarm = (sum(alarm_p) >= length(alarm_p)/2);
end
