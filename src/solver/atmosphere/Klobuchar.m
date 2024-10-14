function [rhos_IF, dIon] = Klobuchar(alpha, beta, rhos, uobs, upvt_LLA0, az, el, eph_dict)
%az 卫星方位角,rad,计算时需要换算为单位pi
%用户地理坐标为经纬高，°，需要换算为单位pi
utime = [uobs.ObsTime];
c = 2.99792458e8;
%计算地心角
PSI = 0.0137./(el/pi + 0.11) - 0.022;
%计算穿刺点地理纬度
PHI = upvt_LLA0(1)/180 + PSI.*cos(az);

PHI(PHI > 0.416)   = 0.416;
PHI(PHI < -0.416) = -0.416;

%计算穿刺点地理经度
lambda = upvt_LLA0(2)/180 + PSI.*sin(az)./cos(PHI*pi);
%计算穿刺点地磁纬度
PHI_m = PHI + 0.064*cos((lambda - 1.617)*pi);
%计算穿刺点地方时
sec = utime;
uweek = floor(sec/(86400*7));
uGPStime = sec - uweek * 86400 * 7;

t = 4.32e4*lambda + uGPStime;

t = t - floor(t/86400)*86400;

%计算投影系数
F = 1 + 16 * (0.53 - el/pi).^3;
%计算电离层延时
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
dIon(1 : length(F)) = F .*(5e-9 + A .* (1 - x.^2/2 + x.^4/24).*(abs(x) < 1.57));
rhos_IF = rhos - c * dIon(1 : length(F));
dIon = dIon';