function [pu, vu, dtu, ddtu, rhor, drhor, H] = lse4pnt(rho, drho, ps, vs, dts, sys)
% estimate user position/velocity/clock bias & drift with pseudoranges
% args  :   1xM     double  rho     [m], bias-fixed pseudorange observables
%           1xM     double  drho    [m/s], bias-fixed pseudorange-rate observables
%           3xM     double  ps      [m], satellite position in ECEF
%           3xM     double  vs      [m/s], satellite velocity in ECEF
%           1xM     double  dts     [s], satellite clock bias
%           1xM     char    sys     satellite system of each obserable
% return:   1x3     double  pu      [m], user position in ECEF
%           1x3     double  vu      [m/s], user velocity in ECEF
%           1xnsys  double  dtu   [s], user clock bias to all N systems
%           double          ddtu    [s/s], use clock drift
%           1xM     double  rhor    [m], pseudorange residual
%           1xM     double  drhor   [m/s], pseudorage rate residual
%           Mx4     double  H       the positioning matrix

    sys(sys=='J') = 'G';            % QZSS is equivalent to GPS in solution
    usys = unique(sys);             % systems involved
    nsys = length(usys);            % number of systems involved
    N = 3 + nsys;                   % 3(position) + nsys(system-time-bias)
    M = length(rho);                % #Satellite used in PVT solution
    if(M < N)
        warning("Unsufficient observables, M=%d", M);
        pu = zeros(1, 3) + NaN;
        vu = zeros(1, 3) + NaN;
        dtu = zeros(1, nsys) + NaN;  ddtu = NaN;
        rhor = NaN; drhor = NaN;
        H = zeros(M, 4) + NaN;
        return;
    end

    pu_k = [0 0 0];                 % Initialized position (ECEF)
    dtu_k = zeros(1, nsys);         % Initialized clock bias
    error_k = 1e5;                  % Position correction 
    iter_max = 100;                 % Iteration allowance
    iter_cnt = 0;
    c = 299792458;                  % Speed of light in GPS
    OMGE = 7.2921151467e-5;

    % Geometry matrix
    H = zeros(M, N);
    b = zeros(M, 1);
    dtu_idx = zeros(1, nsys);

    while error_k>1e-5 && iter_cnt<iter_max
        for i = 1:M
            ps_i = ps(i, :);
            dtu_idx(i) = find(usys==sys(i), 1);
            H(i, 1:3) = -(ps_i-pu_k)/vecnorm(ps_i-pu_k);
            H(i, 3+dtu_idx(i)) = 1;
            r = vecnorm(ps_i-pu_k) + OMGE*(ps_i(1)*pu_k(2)-ps_i(2)*pu_k(1))/c;
            assert(dtu_idx(i) <= length(dtu_k), sprintf("sys=%s, dtu_idx=%d, nsys=%d", sys, dtu_idx, length(dtu_k)));
            b(i) = (r + dtu_k(dtu_idx(i)) - c*dts(i)) - rho(i);
        end
        dx = -H\b;
        error_k = max(abs(dx));
        pu_k = pu_k + dx(1:3)';
        dtu_k = dtu_k + dx(4:end)';
        iter_cnt = iter_cnt+1;
    end
    if(iter_cnt == iter_max)
        warning("PVT solver quits because of #iteration limit.");
    end
    rhor = b;
    pu = pu_k;
    dtu  = dtu_k/c;
    
    Gk = ones(M, 4);
    b = zeros(M, 1);
    for i = 1:M
        ps_i = ps(i, :);
        vi = vs(i, :);
        dp = ps_i - pu;
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