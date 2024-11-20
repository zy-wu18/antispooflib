function p_lla = ecef2lla(p_ecef)

    persistent wgs;
    if(isempty(wgs))
        wgs = wgs84Ellipsoid;
    end
    a = wgs.SemimajorAxis;
    b = wgs.SemiminorAxis;
    e2 = 1-b^2/a^2;
    ep2 = a^2/b^2 - 1;
    
    p = vecnorm(p_ecef(1:2));
    tanu = (p_ecef(3)/p)*(a/b);
    for i = 1:3
        cos2u = 1/(1+tanu^2);
        sin2u = 1 - cos2u;
        tanphi = (p_ecef(3)+ep2*b*sin2u^(3/2))/(p-e2*a*cos2u^(3/2));
        tanu = b/a*tanphi;
    end
    phi = atan(tanphi)/pi*180;
    if(p_ecef(1) >= 0)
        lambda = atan(p_ecef(2)/p_ecef(1))/pi*180;
    elseif(p_ecef(1) < 0 && p_ecef(2) >= 0)
        lambda = atan(p_ecef(2)/p_ecef(1))/pi*180 + 180;
    else
        lambda = atan(p_ecef(2)/p_ecef(1))/pi*180 - 180;
    end
    N = a/sqrt(1-e2*sind(phi)^2);
    if(phi < 90)
        h = p/cosd(phi) - N;
    else
        h = p_ecef(3)/sind(phi) - N + e2*N;
    end
    p_lla = [phi, lambda, h];
end

