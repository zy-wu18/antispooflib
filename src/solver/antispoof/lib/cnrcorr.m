function sdres = cnrcorr(cnrs, keys, cfg)
% CNR correlation detector
% args  :   1xM double  cnrs    [dBHz], carrier-to-noise ratio
%           draimcfg_t  cfg     (optional)customized configuration
% return:   sdres_t     sdres   name, alarm, test statistics and threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(nargin < 2 || isempty(cfg))
        cfg = struct();
        cfg.NumAccum=10;
        cfg.NumPoll=1;
        cfg.Step=10;
        cfg.Thresh=0.5;
    end
    sdres.name = 'CNR Correlation';
    N = cfg.NumAccum;   assert(isscalar(N) && N >= 3);
    K = cfg.NumPoll;    assert(isscalar(K) && K >= 1);
    
    %% Observable accumulation
    persistent i; % function counter
    persistent ig; % global counter
    persistent cnr_arr;
    persistent key_arr;
    
    if(isempty(i) || N~=length(cnr_arr))
        i = 0;
        ig = 0;
        cnr_arr = cell(1, N);
        key_arr = cell(1, N);
    end
    cnr_arr(2:end) = cnr_arr(1:end-1);
    cnr_arr{i+1} = cnrs;
    key_arr(2:end) = key_arr(1:end-1);
    key_arr{i+1} = keys;
    i = mod(i+1, N);
    ig = ig + 1;
    
    sdres.alarm = nan; sdres.T = nan; sdres.gamma = nan;
    if(ig < N || i<cfg.Step-1)
        return;
    end

    %% CNR correlation core
    ikeys = key_arr{1};
    for n = 2:N
        if(length(key_arr{n}) < 4)
            sdres.alarm = 0; sdres.T = nan; sdres.gamma = nan;
            return;
        end
        ikeys = intersect(ikeys, key_arr{n});
    end
    if(length(ikeys) < 2)
        sdres.alarm = nan; sdres.T = nan; sdres.gamma = nan;
        return;
    end

    cnr_vec = zeros(length(ikeys), N);
    for n = 1:N
        cnr = cnr_arr{n};
        for m = 1:length(cnr)
            idx = find(strcmp(ikeys, key_arr{n}(m)), 1);
            if(isempty(idx))
                continue;
            end
            cnr_vec(idx, n) = cnr_arr{n}(m) + 1e-3*randn();
        end
    end
    U = triu(corr(cnr_vec'), 1);
    cross_corr = U(U ~= 0);
    sdres.T = abs(mean(cross_corr));
    sdres.gamma = cfg.Thresh;

    %% Output filtering
    persistent alarm_p;
    if(K ~= length(alarm_p))
        alarm_p = zeros(1, K);
    end
    alarm_t = any(sdres.T > sdres.gamma);
    alarm_p = [alarm_t, alarm_p(1:end-1)];
    sdres.alarm = (sum(alarm_p) >= length(alarm_p)/2);
end