function [el, az] = satelaz_(pu, ps)
% compute satellite azimuth/elevation angle 
% args  :   1x3 double  pu       [m] user position ECEF
%           1x3 double  ps       [m] satellite position ECEF                          
% return:   double      az       [rad] azimuth (0.0<=az<2*pi)
%           double      el       [rad] elevation (-pi/2<=el<=pi/2)
    dxyz = (ps' - pu');
    pu_lla = ecef2lla_(pu);
    lambda = pu_lla(2);
    phi = pu_lla(1);
    denu = [-sind(lambda), cosd(lambda), 0;...
        -sind(phi)*cosd(lambda), -sind(phi)*sind(lambda), cosd(phi); ...
        cosd(phi)*cosd(lambda), cosd(phi)*sind(lambda), sind(phi)]*dxyz;
    az = atan2(denu(1), denu(2));
    el = atan(denu(3)/vecnorm(denu(1:2)));
    
    if(az < 0)
        az = az + 2*pi;
    end
    if(el < 0)
        warning("elv abnormal: elv=%.1f, enu=[%.1e,%.1e,%.1e]", ...
            el, denu(1), denu(2), denu(3));
    end
end