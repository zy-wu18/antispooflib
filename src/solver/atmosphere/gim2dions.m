function [rhos_corr, dions] = gim2dions(ionpath, uobs, upvt)
    persistent GIM;
    persistent t;
    persistent H;
    persistent lats;
    persistent lons;
    persistent vtecs;

    %% Calculate Pierce Lat(phi)/Lon(lambda)
    pu_lla = upvt.PosLLA;
    el = [uobs.El];
    az = [uobs.Az];
    % Calculate the geocentric angle
    psi = 0.0137./(el/pi + 0.11) - 0.022;
    % Calculate the geographic latitude of the puncture point, in semicycle
    phi = pu_lla(1)/180 + psi.*cos(az);
    phi(phi > 0.416) = 0.416;
    phi(phi < -0.416) = -0.416;
    % Calculate the geographic longitude of the puncture point, in semicycle
    lambda = pu_lla(2)/180 + psi.*sin(az)./cos(phi*pi);
    
    %% Load GIM
    if(isempty(GIM))
        GIM = load(ionpath);
        t = GIM.times;
        vtecs = GIM.vtecs;
        lons = GIM.lons;
        lats = GIM.lats;
        H = GIM.H;
    end

    t_tar1 = upvt.Time; % target time
    t_tar1_ymdh = [t_tar1.Year, t_tar1.Month, t_tar1.Day, t_tar1.Hour];
    tidx1 = all(t(:,1:4) == t_tar1_ymdh, 2);
    vtecmap1 = squeeze(vtecs(tidx1, :, :));
    k1 = (60-t_tar1.Minute)/60;

    t_tar2 = upvt.Time + 1/24; % target time
    t_tar2_ymdh = [t_tar2.Year, t_tar2.Month, t_tar2.Day, t_tar2.Hour];
    tidx2 = all(t(:,1:4) == t_tar2_ymdh, 2);
    vtecmap2 = squeeze(vtecs(tidx2, :, :));

    vtecmap = vtecmap1*k1 + vtecmap2*(1-k1);
    
    %% Get STEC satellite by satellite
    M = length(az);
    vtec = zeros(1, M);
    stec = zeros(1, M);
    dions = zeros(1, M);
    meter2tecu = (1.57542e9)^2/(40.3*1e16); % VTEC Map is in L1-meter
    rhos_corr = [uobs.Rho];

    for j = 1:M
        vtec(j) = interp2(lons, lats, vtecmap, lambda(j)*180, phi(j)*180)*meter2tecu;
        %vtec(j) = interp2(lons, lats, vtecmap, pu_lla(2), pu_lla(1))*meter2tecu;
        stec(j) = vtec(j) * 1/sqrt(1 - (6371/(6371+H)*cos(el(j)))^2);
        tecu2meter = 40.3*1e16/uobs(j).Fc^2;
        rhos_corr(j) = uobs(j).Rho - tecu2meter*stec(j);
        dions(j) = tecu2meter*stec(j)/2.99792458e8;
    end
end

