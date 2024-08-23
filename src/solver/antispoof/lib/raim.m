function sdres = raim(rho, ps, cfg)
% RAIM: Receiver autonomous intergrity monitoring with pseudo-range residue
% args  :   1xM double  rho     [m], bias-fixed pseudorange observables
%           3xM double  ps      [m], satellite position in ECEF
%           raimcfg_t   cfg     (optional)customized configuration
% return:   sdres_t     sdres   name, alarm, test statistics and threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(nargin < 3 || isempty(cfg))
        cfg = struct();
        cfg.NumAccum= 1;
        cfg.NumPoll = 1;
        cfg.Pfa = 1e-3;
        cfg.Sigma2 = 300;
    end
    sdres.name = 'RAIM';
    N = cfg.NumAccum;   assert(isscalar(N) && N >= 1);
    K = cfg.NumPoll;    assert(isscalar(K) && K >= 1);

    %% Observable accumulation
    [~, ~, ~, ~, rhor, ~, ~] = lse4pnt(rho, rho+nan, ps, ps+nan);
    persistent i; % thread counter
    persistent M_arr;
    persistent rhor_arr;
    if(isempty(i))
        i = 0;
        M_arr = zeros(1, N);
        rhor_arr = cell(1, N);
    end
    M_arr(i+1) = length(rhor);
    rhor_arr{i+1} = rhor;
    i = mod(i+1, N);

    if(mod(i, N)~= 0 || min(M_arr) <= 4)
        sdres.alarm = nan; sdres.T = nan; sdres.gamma = nan;
        return;
    end

    %% RAIM (Sturza, 1998)
    M = min(M_arr);
    D = zeros(1, N);
    for n = 1:N
        f = rhor_arr{n}(1:M);
        D(n) = f'*f;
    end
    sigma2 = cfg.Sigma2;
    sdres.T = sum(D);
    sdres.gamma = sigma2*chi2inv(1-cfg.Pfa, N*(M-4));
    
    %% Output filtering
    persistent alarm_p;
    if(K ~= length(alarm_p))
        alarm_p = zeros(1, K);
    end
    alarm_t = any(sdres.T > sdres.gamma);
    alarm_p = [alarm_t, alarm_p(1:end-1)];
    sdres.alarm = (sum(alarm_p) >= length(alarm_p)/2);
end
