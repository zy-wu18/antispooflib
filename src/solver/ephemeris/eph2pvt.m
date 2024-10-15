function [ps, vs, dts] = eph2pvt(tsv, eph, iono_opt)
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
        [ps, dts] = neph2pt(tsv, eph, iono_opt);
        dt = 1e-3;
        [ps1, ~ ] = neph2pt(tsv+dt, eph, iono_opt);
        vs = (ps1 - ps)/dt;

    elseif any(eph.sys == ['R' 'S'])
        [ps, vs, dts] = geph2pvt(tsv, eph);
    else
        error("InvalidObservableType");
    end
end

