function [p, dp] = gephiter(p0, dp0, ddp, tspan)
% GEPHITER GLONSS ephemeris iterator
% args  :   3x1 double  p0      [m], satellite ECEF position at t0
%           3x1 double  dp0     [m/s], satellite ECEF velocity at t0
%           3x1 double  ddp0    [m/s^2], satellite ECEF acceleration at t0
%           double      tspan   [s], iteration time span
% return:   3x1 double  p       [m], satellite ECEF position at (t0+tspan)
%           3x1 double  dp      [m/s], satellite ECEF velocity at (t0+tspan)
% notes :   based on ode45 non-stiff differential equatio solver

    [mu, omega, ~] = gnssconst('R');
    a_e = 6378136; % earth long axis, [m]
    J2 = 1082625.7e-9; % Gravity field 2nd resonence coeff
    alpha = 3/2*J2*mu*(a_e^2);
    
    f = @(t, y)[y(4); y(5); y(6);
        -mu*y(1)/norm(y(1:3))^3 + alpha/(norm(y(1:3))^5)*y(1)*(1-5*y(3)^2/norm(y(1:3))^2) + ddp(1) + omega^2*y(1) + 2*omega*y(5);
        -mu*y(2)/norm(y(1:3))^3 + alpha/(norm(y(1:3))^5)*y(2)*(1-5*y(3)^2/norm(y(1:3))^2) + ddp(2) + omega^2*y(2) - 2*omega*y(4);
        -mu*y(3)/norm(y(1:3))^3 + alpha/(norm(y(1:3))^5)*y(3)*(3-5*y(3)^2/norm(y(1:3))^2) + ddp(3)];
    y0 = [p0; dp0];
    opt = odeset("RelTol", 1e-2, "AbsTol", 1e-4);
    [t, y] = ode45(f, [0, tspan], y0, opt);
    p = y(end, 1:3);
    dp = y(end, 4:6);
end

