function dt = eph2clk(tsv, eph)
%EPH2CLK GPS/BDS/Galileo eph to satellite clock fix
% args  :   double  tsv     [s], transmitting time
%           eph_t   eph     ephemeris struct.
% return:   double  dt      [s], satellite clock fix
    sys = eph.sys;
    if(any(sys == 'GJEC'))
        t = tsv - eph.Toc_TOW;
        ts = tsv - eph.Toc_TOW;
        for i = 0:2
            t = ts - (eph.af0 + eph.af1*t + eph.af2*t*t);
        end
        dt = eph.af0 + eph.af1*t + eph.af2*t*t;
    elseif(any(sys == 'RS'))
        t = tsv - eph.toe;
        ts = tsv - eph.toe;
        for i = 0:2
            t = ts - (-eph.TauN + eph.GammaN*t);
        end
        dt = -eph.TauN + eph.GammaN;
    else
        warning('eph2clk: unknown system %s', sys);
        dt = nan;
    end
end

