function [pu, vu, dtu, ddtu, rhor, drhor, H] = lse7pnt(drho, ps, vs)
% estimate user position/velocity & clock drift with pseudorange rates
% args  :   Mx1 double  drho    [m/s], bias-fixed pseudorange-rate observables
%           Mx3 double  ps      [m], satellite position in ECEF
%           Mx3 double  vs      [m/s], satellite velocity in ECEF
% return:   1x3 double  pu      [m], user position in ECEF
%           1x3 double  vu      [m/s], user velocity in ECEF
%           double      dtu     (reserved)[s], always 0
%           double      ddtu    [s/s], use clock drift
%           Mx1 double  rhor    [m], pseudorange residual
%           Mx1 double  drhor   [m/s], pseudorage rate residual
%           Mx7 double  H       the positioning matrix

    M = length(drho);                % #Satellite used in PVT solution
    assert(isnumeric(drho) & iscolumn(drho));
    assert(isnumeric(ps) & size(ps, 2) == 3 & size(ps, 1) == M);
    assert(isnumeric(vs) & size(vs, 2) == 3 & size(vs, 1) == M);
    
    if(M < 7)
        warning("Unsufficient observables, M=%d", M);
        pu = nan(1, 3);
        vu = nan(1, 3);
        dtu = [];
        ddtu = nan;
        rhor = nan(M, 1);
        drhor = nan(M, 1);
        H = nan(M, 7);
        return;
    end

    pu_k = [0 0 0];                 % Initialized position (ECEF), [m]
    vu_k = [0 0 0];                 % Initialized velocity (ECEF), [m]
    ddtu_k = 0;                     % Initialized clock drift, [s/s], scale=c
    error_k = 1e5;                  % Position correction 
    iter_max = 100;                 % Iteration allowance
    iter_cnt = 0;
    c = 299792458;                  % Speed of light in GPS
    
    % Geometry matrix
    H = zeros(M, 7);
    b = zeros(M, 1);
    
    while error_k>1e-5 && iter_cnt<iter_max
        for i = 1:M
            pi = ps(i, :);
            vi = vs(i, :);
            dp_norm = (pi - pu_k) / vecnorm(pi - pu_k);
            dv_norm = (vi - vu_k) / vecnorm(pi - pu_k);
            H(i, 1:3) = cross(dp_norm, cross(dp_norm, dv_norm));
            H(i, 4:6) = -1.0*dp_norm;
            H(i, 7) = 1.0;
            b(i) = (dot(vi - vu_k, dp_norm) + ddtu_k) - drho(i);
        end
        dx = -H\b;
        error_k = max(abs(dx));
        pu_k = pu_k + dx(1:3)';
        vu_k = vu_k + dx(4:6)';
        ddtu_k = ddtu_k + dx(7);
        iter_cnt = iter_cnt+1;
    end
    if(iter_cnt == iter_max)
        warning("PVT solver quits because of #iteration limit.");
    end
    
    pu = pu_k;
    vu = vu_k;
    dtu = [];
    ddtu = ddtu_k/c;
    rhor = zeros(1, M) + NaN;
    drhor = b;
end