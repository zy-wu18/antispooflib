function [ps, vs, dts] = geph2pvt(tsv, eph)
% Calculate the position, velocity and clock bias of a GLONASS satellite
% args  :   double      tsv     transmit time, tsv = tlatch - rho/c
%           eph_t       eph     GLONASS ephemeris data struct
% return:   1x3 double  ps      [m], satellite ECEF position [x, y, z]
%           double      vs      [m/s], satellite ECEF velocity
%           double      dts     [s], satellite clock fix
% notes :   Transfer GLONASS ephemeris at transmt time to satellite clock 
%           fix <dts>, satellite position <ps> and velocity <vs>

    t = tsv - eph.toe;
    dts = -eph.TauN + eph.GammaN*t;

    pos = 1e3*[eph.X; eph.Y; eph.Z];
    vel = 1e3*[eph.Xvel; eph.Yvel; eph.Zvel];
    acc = 1e3*[eph.Xacc; eph.Yacc; eph.Zacc]; % Gravity field 2nd resonance
    iterstep = (2*(t>=0)-1)*60;
    while(abs(t) > 1e-9)
        if(abs(t) < 60)
            iterstep = t;
        end
        [pos, vel] = glorbit(iterstep, [pos; vel], acc);
        t = t - iterstep;
    end
    ps = pos;
    vs = vel;
end

function [ps, vs] = glorbit(t, x, acc)
    [mu, omega] = gnssconst('R');
    a_e = 6378136; % earth long axis, [m]
    J2 = 1082625.7e-9; % Gravity field 2nd resonence coeff
    alpha = 3/2*J2*mu*(a_e^2);
    f = @(y)[y(4); y(5); y(6);
        -mu*y(1)/norm(y(1:3))^3 + alpha/(norm(y(1:3))^5)*y(1)*(1-5*y(3)^2/norm(y(1:3))^2) + acc(1) + omega^2*y(1) + 2*omega*y(5);
        -mu*y(2)/norm(y(1:3))^3 + alpha/(norm(y(1:3))^5)*y(2)*(1-5*y(3)^2/norm(y(1:3))^2) + acc(2) + omega^2*y(2) - 2*omega*y(4);
        -mu*y(3)/norm(y(1:3))^3 + alpha/(norm(y(1:3))^5)*y(3)*(3-5*y(3)^2/norm(y(1:3))^2) + acc(3)];
    
    k = zeros(6, 4);
    k(:,1) = f(x); w = x + k(:,1)*t/2;
    k(:,2) = f(w); w = x + k(:,2)*t/2;
    k(:,3) = f(w); w = x + k(:,3)*t;
    k(:,4) = f(w);
    weight = [1/6; 2/6; 2/6; 1/6]*t;
    dx = k*weight;
    ps = x(1:3) + dx(1:3);
    vs = x(4:6) + dx(4:6);
end
