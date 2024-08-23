function [ps, vs, dts] = eph2pvt(tsv, eph)
% Calculate the position, velocity and clock bias of a single satellite
% args  :   double      tsv     transmit time, tsv = tlatch - rho/c
%           eph_t       eph     ephemeris data struct
% return:   1x3 double  ps      [m], satellite ECEF position [x, y, z]
%           1x3 double  vs      [m/s], satellite ECEF velocity [x, y, z]
%           double      dts     [s], satellite clock fix
% notes :   GPS, GAL, QZSS, BDS = calculate satellite position, then
%           calculate velocity as (pos@(t+dt)-pos@(t))/dt, dt=1ms
%           GLONASS, SBAS = calculate both position and velocity using
%           runge_katta_4 iterator.

    if any(eph.sys == ['G' 'E' 'J' 'C'])
        [ps, dts] = neph2pos(tsv, eph);
        dt = 1e-3;
        [ps1, ~ ] = neph2pos(tsv+dt, eph);
        vs = (ps1 - ps)/dt;

    elseif any(eph.sys == ['R' 'S'])
        % Satellite clock bias fix
        t_toc = tsv - eph.toe; % t \approx tsv
        dtsv = -eph.TauN + eph.GammaN*t_toc; % Satellite clock-bias correction
        t = tsv - dtsv; % satellite clock fixed
        tk = t - eph.toe; % Time elapsed since ephemeris updated at Toe
        
        %% User Runge-Kutta (4th order) to solve GLONASS satellite PVT
        p0  = 1e3*[eph.X; eph.Y; eph.Z];
        dp0 = 1e3*[eph.Xvel; eph.Yvel; eph.Zvel];
        ddp = 1e3*[eph.Xacc; eph.Yacc; eph.Zacc]; % Gravity field 2nd resonance
        [ps, vs] = gephiter(p0, dp0, ddp, tk);
    
        %% Satellite time bias
        dts = dtsv;
    else
        error("InvalidObservableType");
    end
    
end

