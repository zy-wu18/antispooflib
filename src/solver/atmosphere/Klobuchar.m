function [rhos_K, dIon] = klobuchar(rhos, uobs, upvt_LLA0, az, el, eph_dict)
% Calculate ionospheric delay using Klobuchar model
% args  :   4x1 double  alpha     ionospheric parameters in broadcast ephemeris
%           4x1 double  beta      ionospheric parameters in broadcast ephemeris
%           1xM double  rhos      [m], satellite pseodoranges
%           1xM struct  uobs      satellites observation used in PVT calculation
%           1x1 struct  upvt      previous receiver PVT observation
%           1x3 double  upvt_LLA0 [°,°,m] receiver LLA position
%           1xM double  az,el     [rad], satellite direction angle and elevation angles
%           dictionary  eph_dict  ephemeris data struct
% return:   1xM double  dIon     [s], ionospheric delay
%           1xM double  rhos_IF   [m], satellite pseodoranges corrected
% Notes:    all angles used should be converted to semi-circle(pi) for calculation
%           TGD has been included in rhos_IF

M = length(uobs);
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
utime = zeros(1, M);
for k=1:M
    [uGPStime, ~] = Utc2Gps([uobs(k).Time]);
    if(strcmp(uobs(k).Sys, 'G'))
        utime(k) = uGPStime(2);
    elseif(strcmp(uobs(k).Sys, 'C'))
        utime(k) = uGPStime(2) - 14;
    end
end

t = 4.32e4*lambda + utime;
t = t - floor(t/86400)*86400;
% Calculate projection coefficient
F = 1 + 16 * (0.53 - el/pi).^3;
% Calculate ionospheric delay
B = zeros(1, M);
A = zeros(1, M);
x = zeros(1, M);

for m = 1:M
    alpha = eph_dict(['iono_', uobs(1).Sys]).alpha;
    beta = eph_dict(['iono_', uobs(1).Sys]).beta;
    B(m) = (PHI_m(m).^(0:3)) * beta;
    if(B(m) < 72000)
        B(m) = 72000;
    end
    x(m) = 2*pi*(t(m) - 50400)./B(m);
    
    A(m) = (PHI_m(m).^(0:3)) * alpha;
    if(A(m) < 0)
        A(m) = 0;
    end
end
dIon = F .*(5e-9 + A .* (1 - x.^2/2 + x.^4/24).*(abs(x) < 1.57));
rhos_K = rhos - c*dIon;

% TGD correction
for k = 1:M
    key = sprintf("%c%02d", uobs(k).Sys, uobs(k).PRN);
    if(eph_dict.isKey(key) && ~isempty(eph_dict(key).TGD))
        TGD = eph_dict(key).TGD;
        rhos_K(k) = rhos_K(k) - c*TGD(1);
    end
end 