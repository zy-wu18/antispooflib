function sdres = draim2nd(drho, ps, vs, cfg)
% Doppler RAIM detector 2nd model
% args  :   1xM double  drho    [m/s], bias-fixed pseudorange-rate
%           3xM double  ps      [m], satellite position in ECEF
%           3xM double  vs      [m/s], satellite velocity in ECEF
%           draimcfg_t  cfg     (optional)customized configuration
% return:   sdres_t     sdres   name, alarm, test statistics and threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(nargin < 3 || isempty(cfg))
        cfg = struct();
        cfg.NumAccum= 4;
        cfg.NumPoll = 5;
        cfg.Pfa = 1e-3;
    end
    sdres.name = 'DRAIM2nd';
    N = cfg.NumAccum;   assert(isscalar(N) && N >= 3);
    K = cfg.NumPoll;    assert(isscalar(K) && K >= 1);

    %% Observable accumulation
    [~, ~, ~, ~, ~, drhor, H7] = lse7pnt(drho, ps, vs);
    persistent i; % thread counter
    persistent M_arr;
    persistent H7_arr;
    persistent drhor_arr;

    if(isempty(i))
        i = 0;
        M_arr = zeros(1, N);
        H7_arr = cell(1, N);
        drhor_arr = cell(1, N);
    end
    M_arr(i+1) = length(drhor);
    H7_arr{i+1} = H7;
    drhor_arr{i+1} = drhor;
    i = mod(i+1, N);

    if(mod(i, N)~= 0 || min(M_arr) <= 7)
        sdres.alarm = nan; sdres.T = nan; sdres.gamma = nan;
        return;
    end
    
    %% GLRT Doppler RAIM Core (Zhou, et.al, 2023)
    M = min(M_arr);
    UT = zeros(N*(M-7), N*M);
    GT = zeros(3*(M-7), N*(M-7));
    r = zeros(N*M, 1);
    for n = 1:N
        Hn = H7_arr{n}(1:M, :); % truncate extra hn in window 
        Sn = eye(M) - Hn/(Hn'*Hn)*Hn';
        [Un, ~, ~] = svd(Sn);
        UT(((n-1)*(M-7)+1):n*(M-7), (n-1)*M+1:n*M) = Un(:, 1:M-7)';
        GT(1:M-7, ((n-1)*(M-7)+1):n*(M-7))             = eye(M-7);
        GT(1*(M-7)+1:2*(M-7), ((n-1)*(M-7)+1):n*(M-7)) = (n-1)*eye(M-7);
        GT(2*(M-7)+1:3*(M-7), ((n-1)*(M-7)+1):n*(M-7)) = (n-1)^2*eye(M-7);
        r(((n-1)*M+1):n*M) = drhor_arr{n}(1:M);
    end
    G = GT';
    alpha = (GT*G)\GT*UT*r;
    sigma2_f0 = (UT*r)'*(UT*r);
    sigma2_f1 = (UT*r-G*alpha)'*(UT*r-G*alpha);
    sdres.T = (N-3)/3*(sigma2_f0 - sigma2_f1)/sigma2_f1;
    sdres.gamma = finv(1-cfg.Pfa, 3*(M-7), (N-3)*(M-7));

    %% Output filtering
    persistent alarm_p;
    if(K ~= length(alarm_p))
        alarm_p = zeros(1, K);
    end
    alarm_t = any(sdres.T > sdres.gamma);
    alarm_p = [alarm_t, alarm_p(1:end-1)];
    sdres.alarm = (sum(alarm_p) >= length(alarm_p)/2);
end

