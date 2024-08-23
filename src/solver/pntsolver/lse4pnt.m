function [pu, vu, dtu, ddtu, rhor, drhor, H] = lse4pnt(rho, drho, ps, vs)
% estimate user position/velocity/clock bias & drift with pseudoranges
% args  :   1xM double  rho     [m], bias-fixed pseudorange observables
%           1xM double  drho    [m/s], bias-fixed pseudorange-rate observables
%           3xM double  ps      [m], satellite position in ECEF
%           3xM double  vs      [m/s], satellite velocity in ECEF
% return:   1x3 double  pu      [m], user position in ECEF
%           1x3 double  vu      [m/s], user velocity in ECEF
%           double      dtu     [s], user clock bias
%           double      ddtu    [s/s], use clock drift
%           1xM double  rhor    [m], pseudorange residual
%           1xM double  drhor   [m/s], pseudorage rate residual
%           Mx4 double  H       the positioning matrix

    M = length(rho);                % #Satellite used in PVT solution
    if(M < 4)
        warning("Unsufficient observables, M=%d", M);
        pu = zeros(1, 3) + NaN;
        vu = zeros(1, 3) + NaN;
        dtu = NaN;  ddtu = NaN;
        rhor = NaN; drhor = NaN;
        H = zeros(M, 4) + NaN;
        return;
    end

    p_k = [0 0 0];                  % Initialized position (ECEF)
    dt_k = 0;                       % Initialized clock bias
    error_k = 1e5;                  % Position correction 
    iter_max = 100;                 % Iteration allowance
    iter_cnt = 0;
    c = 299792458;                  % Speed of light in GPS
    
    % Geometry matrix
    H = zeros(M, 4);
    b = zeros(M, 1);
    
    while error_k>1e-5 && iter_cnt<iter_max
        for i = 1:M
            pi = ps(i, :);
            dp = pi - p_k;
            dp_norm = vecnorm(dp);
            H(i, 1:3) = -dp/dp_norm;
            H(i, 4) = 1;
            b(i) = dp_norm + dt_k - rho(i);
        end
        dx = -H\b;
        error_k = max(abs(dx));
        p_k = p_k + dx(1:3)';
        dt_k = dt_k + dx(4);
        iter_cnt = iter_cnt+1;
    end
    if(iter_cnt == iter_max)
        warning("PVT solver quits because of #iteration limit.");
    end
    rhor = b;
    pu = p_k;
    dtu  = dt_k/c;
    
    Gk = ones(M, 4);
    b = zeros(M, 1);
    for i = 1:M
        pi = ps(i, :);
        vi = vs(i, :);
        dp = pi - pu;
        dp_norm = norm(dp,'fro');
        los = dp./dp_norm; % Light of sight vector
        Gk(i, 1:3) = -los;
        b(i, 1) = -dot(vi, los) + drho(i);
    end
    drhor = b;
    dx = Gk\b;
    vu = dx(1:3, 1)';
    ddtu = dx(4, 1);
end