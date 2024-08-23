function [pu, vu, dtu, ddtu, rhor, drhor, H] = lse7pnt(drho, ps, vs)
% estimate user position/velocity & clock drift with pseudorange rates
% args  :   1xM double  drho    [m/s], bias-fixed pseudorange-rate observables
%           3xM double  ps      [m], satellite position in ECEF
%           3xM double  vs      [m/s], satellite velocity in ECEF
% return:   1x3 double  pu      [m], user position in ECEF
%           1x3 double  vu      [m/s], user velocity in ECEF
%           double      dtu     (reserved)[s], always 0
%           double      ddtu    [s/s], use clock drift
%           1xM double  rhor    [m], pseudorange residual
%           1xM double  drhor   [m/s], pseudorage rate residual
%           Mx7 double  H       the positioning matrix

    M = length(drho);                % #Satellite used in PVT solution
    if(M < 7)
        warning("Unsufficient observables, M=%d", M);
        pu = zeros(1, 3) + NaN;
        vu = zeros(1, 3) + NaN;
        dtu = 0.0;
        ddtu = NaN;
        rhor = zeros(1, M) + NaN;
        drhor = zeros(1, M) + NaN;
        H = zeros(M, 7) + NaN;
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
    dtu = 0.0;
    ddtu = ddtu_k/c;
    rhor = zeros(1, M) + NaN;
    drhor = b;
end