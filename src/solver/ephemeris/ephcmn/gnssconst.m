function [mu, o_dot_e, F] = gnssconst(sys)
% Get GNSS constants
% args  :   char    sys     I   'C','E','G','J','R'
% return:   double  mu      [m^3/s^2], Universal gravitational constant
%           double  O_dot_e [rad/s], earth angle speed
%           double  F       [sec/(m^(1/2))] Theory of Relativity fix
    switch(sys)
        case 'C'
            mu      = 3.986004418e14;
            o_dot_e = 7.2921150e-5;
        case 'E'
            mu      = 3.986004418e14; 
            o_dot_e = 7.2921151467e-5; 
        case {'G', 'J'}
            mu      = 3.986005e14;    
            o_dot_e = 7.2921151467e-5;
        case 'R'
            mu      = 3.9860044e14;
            o_dot_e = 7.2921150e-5;
        otherwise
            mu      = 3.986005e14;    
            o_dot_e = 7.2921151467e-5;
    end
    c = 2.99792458e8; % [m/s], speed of light
    F = -2*sqrt(mu)/(c^2); % [sec/(m^(1/2))] Theory of Relativity
end

