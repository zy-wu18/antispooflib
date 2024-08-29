function [ps, vs, dts] = geph2pvt(tsv, eph)
% Calculate the position, velocity and clock bias of a GLONASS satellite
% args  :   double      tsv     transmit time, tsv = tlatch - rho/c
%           eph_t       eph     GLONASS ephemeris data struct
% return:   1x3 double  ps      [m], satellite ECEF position [x, y, z]
%           double      vs      [m/s], satellite ECEF velocity
%           double      dts     [s], satellite clock fix
% notes :   Transfer GPS, GAL, QZSS, BDS ephemeris at transmt time to 
%           satellite clock fix <dts> and satellite position <ps>
    tk = tsv - eph.toe;
    p0  = 1e3*[eph.X; eph.Y; eph.Z];
    dp0 = 1e3*[eph.Xvel; eph.Yvel; eph.Zvel];
    ddp = 1e3*[eph.Xacc; eph.Yacc; eph.Zacc]; % Gravity field 2nd resonance
    [ps, vs] = gephiter(p0, dp0, ddp, tk);
    
    dts = eph2clk(tsv, eph);
end

