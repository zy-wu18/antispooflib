function [ps, dts] = neph2pt(tsv, eph, iono_opt)
% Calculate the position, velocity and clock bias of a G/C/J/E satellite
% args  :   double      tsv     transmit time, tsv = tlatch - rho/c
%           eph_t       eph     G/C/J/E ephemeris data struct
% return:   1x3 double  ps      [m], satellite ECEF position [x, y, z]
%           double      dts     [s], satellite clock fix
% notes :   Transfer GPS, GAL, QZSS, BDS ephemeris at transmt time to 
%           satellite clock fix <dts> and satellite position <ps>

    assert(any(eph.sys == ['G' 'E' 'J' 'C']));
    [mu, Omega_dot_ie, F] = gnssconst(eph.sys);

    %% Load ephemeris constants
    toc  = eph.Toc_TOW; af0 = eph.af0;  af1     = eph.af1;      af2     = eph.af2;
    IODE = eph.IODE;    Crs = eph.Crs;  Delta_n = eph.Delta_n;  M_0     = eph.M_0;
    Cuc  = eph.Cuc;     e   = eph.e;    Cus     = eph.Cus;      sqrt_a  = eph.sqrt_a;
    Toe  = eph.Toe;     Cic = eph.Cic;  Omega_0 = eph.Omega_0;  Cis     = eph.Cis;
    i_0  = eph.i_0;     Crc = eph.Crc;  omega   = eph.omega;    Omega_dot=eph.Omega_dot;
    I_dot= eph.I_dot;   TGD = eph.TGD;  WN      = eph.WN;       Toes    = eph.Toes;

    %% Begin to calculate satellite position/velocity
    % Kepler
    tk = tsv - Toe; % ephemeris dead reaconing duration
    n = sqrt(mu./sqrt_a.^6) + Delta_n; % Average velocity
    Mk = M_0 + n.*tk;
    Ek1 = Mk; Ek = Ek1 - (Ek1 - e*sin(Ek1) - Mk)/(1 - e*cos(Ek1));
    while(abs(Ek-Ek1) > 1e-10)
        Ek1 = Ek; Ek = Ek1 - (Ek1 - e*sin(Ek1) - Mk)/(1 - e*cos(Ek1));
    end
    % True anomaly
    vk = atan2(sqrt(1-e.^2).*sin(Ek), cos(Ek) - e);
    % Ascending pitch angle and its rate of change
    phik = vk + omega;
    % Calculate perturbation correction terms for velocity direction, radial and inclination direction
    delta_uk = Cus.*sin(2*phik) + Cuc.*cos(2*phik); % velocity correction
    delta_rk = Crs.*sin(2*phik) + Crc.*cos(2*phik); % radius correction
    delta_ik = Cis.*sin(2*phik) + Cic.*cos(2*phik); % pitch correction
    % Calculate the perturbation corrected elevation angle, satellite vector diameter, orbit inclination angle, and their rate of change
    uk = phik + delta_uk; % Corrected ascending pitch angle
    rk = (sqrt_a.^2)*(1-e.*cos(Ek)) + delta_rk; % Corrected radius
    ik = i_0 + I_dot.*tk + delta_ik; % Corrected pitch angle
    % Calculate the position and rate of change of satellites in the orbital plane coordinate system
    xk = rk.*cos(uk);
    yk = rk.*sin(uk);
    
    %% Dump pos/vel during the observation of a satellite
    % The 3D position/velocity of satellites in ECEF coordinates
    ps = zeros(length(tsv), 3);
    
    % Corrected ascending intersection longitude
    if eph.sys == 'C' && (eph.PRN <= 5 || eph.PRN >= 59) % inertial coordinate
        Omegak = Omega_0 + Omega_dot.*tk - Omega_dot_ie.*Toes; 
        xg = xk.*cos(Omegak) - yk.*cos(ik).*sin(Omegak);
        yg = xk.*sin(Omegak) + yk.*cos(ik).*cos(Omegak);
        zg = yk.*sin(ik);
        ox = Omega_dot_ie*tk;
        oz = -5/180*pi;
        ps(:, 1) = xg.*cos(ox) + yg*sin(ox)*cos(oz) + zg*sin(ox)*sin(oz);
        ps(:, 2) =-xg.*sin(ox) + yg*cos(ox)*cos(oz) + zg*cos(ox)*sin(oz);
        ps(:, 3) =-yg.*sin(oz) + zg*cos(oz);
    else % earth-fixed coo
        Omegak = Omega_0 + (Omega_dot-Omega_dot_ie).*tk - Omega_dot_ie.*Toes; 
        ps(:, 1) = xk.*cos(Omegak) - yk.*cos(ik).*sin(Omegak);
        ps(:, 2) = xk.*sin(Omegak) + yk.*cos(ik).*cos(Omegak);
        ps(:, 3) = yk.*sin(ik);
    end
    
    % Satellite time
    dtsv = eph2clk(tsv, eph);
    dtsr = F*e*sqrt_a*sin(Ek); % relativity fix
    if (strcmp(iono_opt,'IonoFree') || strcmp(eph.sys, 'C'))
        dts = dtsv + dtsr;
    else
        dts = dtsv + dtsr - TGD(1);
    end
end

