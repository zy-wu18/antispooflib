function [el, az] = satelaz(pu, ps)
% compute satellite azimuth/elevation angle 
% args  :   1x3 double  pu       [m] user position ECEF
%           1x3 double  ps       [m] satellite position ECEF                          
% return:   double      az       [rad] azimuth (0.0<=az<2*pi)
%           double      el       [rad] elevation (-pi/2<=el<=pi/2)
    lla = ecef2lla(pu);
    wgs84 = wgs84Ellipsoid;
    [eus, nus, uus] = ecef2enu(ps(1), ps(2), ps(3), lla(1), lla(2), lla(3), wgs84);
    az = atan2(eus, nus);
    el = atan(uus/sqrt(eus^2 + nus^2));
    if(az < 0)
        az = az + 2*pi;
    end
    if(el < 0)
        warning("elv abnormal: elv=%.1f, enu=[%.1e,%.1e,%.1e]", el, eus, nus, uus);
    end
end

