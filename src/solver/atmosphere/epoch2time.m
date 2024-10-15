function [utime, usec] = epoch2time(u_datetime)
% convert datetime to TOW
% args  :   datetime    u_datetime      query UTC time
% return:   int         utime           [s],integer part of TOW
%           double      usec            [s],fractional part of TOW
doy = [1,32,60,91,121,152,182,213,244,274,305,335];
uyear = year(u_datetime);
umon = month(u_datetime);
uday = day(u_datetime);
    
if (uyear<1970 || uyear>2099 || umon<1 || umon>12)
    utime = 0;
else
    days = (uyear - 1970)*365 + (uyear - 1969)/4 + doy(umon) + uday - 2 +...
        ((mod(uyear, 4)==0) && (umon>=3));
    usec = second(u_datetime);
    utime = days*86400 + hour(u_datetime)*3600 + minute(u_datetime)*60 + floor(usec);
    usec = usec - floor(usec);
end