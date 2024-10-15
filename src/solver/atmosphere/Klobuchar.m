function [rhos_IF, dIon] = klobuchar(alpha, beta, rhos, uobs, upvt_LLA0, az, el, eph_dict)
% Calculate ionospheric delay using Klobuchar model
% args  :   4x1 double  alpha     ionospheric parameters in broadcast ephemeris
%           4x1 double  beta      ionospheric parameters in broadcast ephemeris
%           1xM double  rhos      [m], satellite pseodoranges
%           1xM struct  uobs      satellites observation used in PVT calculation
%           1x1 struct  upvt      previous receiver PVT observation
%           1x3 double  upvt_LLA0 [°,°,m] receiver LLA position
%           1xM double  az,el     [rad], satellite direction angle and elevation angles
%           dictionary  eph_dict  ephemeris data struct
% return:  32x1 double  dIon     [s], ionospheric delay
%           1xM double  rhos_IF   [m], satellite pseodoranges corrected
% Notes:    all angles used should be converted to semi-circle(pi) for calculation


prn_num = length(uobs);
c = 2.99792458e8;
% Calculate the geocentric angle
PSI = 0.0137./(el/pi + 0.11) - 0.022;
% Calculate the geographic latitude of the puncture point
PHI = upvt_LLA0(1)/180 + PSI.*cos(az);
PHI(PHI > 0.416)   = 0.416;
PHI(PHI < -0.416) = -0.416;
% Calculate the geographic longitude of the puncture point
lambda = upvt_LLA0(2)/180 + PSI.*sin(az)./cos(PHI*pi);
% Calculate the geomagnetic latitude of the puncture point
PHI_m = PHI + 0.064*cos((lambda - 1.617)*pi);

% Calculate the puncture point location time
if(strcmp(uobs(1).Sys, 'C'))
    [t0, ~] = epoch2time(datetime(2006, 1, 1, 0, 0, 0));
elseif(strcmp(uobs(1).Sys, 'G'))
    [t0, ~] = epoch2time(datetime(1980, 1, 6, 0, 0, 0));
end
sec = zeros(1, prn_num);
for k = 1:prn_num
    [sec(k), sec_ms] = epoch2time(datetime(uobs(k).Time));
    sec(k) = sec(k) - t0 + sec_ms;
end
uweek = floor(sec/(86400*7));
uGPStime = sec - uweek * 86400 * 7;

t = 4.32e4*lambda + uGPStime;
t = t - floor(t/86400)*86400;
% Calculate projection coefficient
F = 1 + 16 * (0.53 - el/pi).^3;
% Calculate ionospheric delay
B = 0;
for k = 1:4
    B = B + beta(k)*PHI_m.^(k - 1);
end
B(B < 72000) = 72000;
x = 2*pi*(t - 50400)./B;
A = 0;
for k = 1:4
    A = A + alpha(k)*PHI_m.^(k - 1);
end
A(A < 0) = 0;

dIon = zeros(1, 32);
dIon(1 : prn_num) = F .*(5e-9 + A .* (1 - x.^2/2 + x.^4/24).*(abs(x) < 1.57));
rhos_IF = rhos - c * dIon(1 : prn_num);
% TGD correction for BDS B1I
for k = 1:prn_num
    if(strcmp(uobs(k).Sys, 'C'))
        if(uobs(k).PRN<10)
            str = [uobs(k).Sys, '0', num2str(uobs(k).PRN)];
        else
            str = [uobs(k).Sys, num2str(uobs(k).PRN)];
        end
        TGD = eph_dict(str).TGD;
        rhos_IF(k) = rhos_IF(k) - c*TGD(1);
    end
end
dIon = dIon';