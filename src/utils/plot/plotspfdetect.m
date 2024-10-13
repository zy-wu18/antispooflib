function fig_sd = plotspfdetect(sdres_seq, obs_rate)
    nmthd = size(sdres_seq, 1); % number of spoofing detector methods
    alarm = reshape([sdres_seq.alarm], nmthd, []);
    T = reshape([sdres_seq.T], nmthd, []);
    gamma = reshape([sdres_seq.gamma], nmthd, []);
    
    if(isempty(sdres_seq))
        fprintf("plotspfdetect: no spoofing-detecion method to draw.\n");
        return
    end
    fig_h = min(512, 256*nmthd);
    fig_sd = figure('Position', [120, 120, 512, fig_h]);
    L = length(alarm);
    tobs = (1/obs_rate)*(1:L);

    for i = 1:nmthd
        subplot(nmthd, 1, i);
        tcurve = T(i, :)./gamma(i, :);
        vld = ~isnan(tcurve);
        if(~any(vld==true))
            continue;
        end
        semilogy(tobs, tcurve, '-o', 'LineWidth', 1.0);
        hold on;
        title(sdres_seq(i, 1).name);
        semilogy([min(tobs(vld)), max(tobs(vld))], [1 1], 'LineWidth', 1.0);
        vld = ~isnan(alarm(1, :)); % only emphasize the alarmed cases
        if(sum(vld) <= 1)
            continue;
        end
        semilogy(tobs, alarm(i, :), 'r', 'LineWidth', 1.0, 'Marker', 'square', 'MarkerFaceColor','r');
        grid on;
        ymax = max(max(T(i, :)./gamma(i, :)*2), 1.5);
        ymin = min(min(T(i, :)./gamma(i, :)/2), 0.2);
        axis([min(tobs(vld)), max(tobs(vld)), ymin, ymax]);
    end
end
