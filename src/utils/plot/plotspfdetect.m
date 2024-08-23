function fig_sd = plotspfdetect(sdres_seq, t)
    nmthd = size(sdres_seq, 1); % number of spoofing detector methods
    alarm = reshape([sdres_seq.alarm], nmthd, []);
    T = reshape([sdres_seq.T], nmthd, []);
    gamma = reshape([sdres_seq.gamma], nmthd, []);
    
    fig_h = min(512, 256*nmthd);
    fig_sd = figure('Position', [120, 120, 512, fig_h]);
    rec_t_utc_sec = [3600, 60, 1]*[t.Hour; t.Minute; t.Second];
    dt = rec_t_utc_sec - rec_t_utc_sec(1);
    if(isempty(sdres_seq))
        fprintf("plotspfdetect: no spoofing-detecion method to draw.\n");
        return
    elseif(all(isnan(dt)))
        fprintf("plotspfdetect: all t is NaN, no time result.\n");
        return
    end

    for i = 1:nmthd
        subplot(nmthd, 1, i);
        tcurve = T(i, :)./gamma(i, :);
        vld = and(~isnan(tcurve), ~isnan(dt));
        semilogy(dt(vld), tcurve(vld), 'LineWidth', 1.0);
        hold on;
        semilogy([min(dt(vld)), max(dt(vld))], [1 1], 'LineWidth', 1.0);
        vld = (alarm(i, :)==1); % only emphasize the alarmed cases
        semilogy(dt(vld), alarm(i, vld), 'ro', 'LineWidth', 1.0);
        grid on;
        ymax = max(max(T(i, :)./gamma(i, :)*2), 1.5);
        ymin = min(min(T(i, :)./gamma(i, :)/2), 0.2);
        axis([min(dt(vld)), max(dt(vld)), ymin, ymax]);
        title(sdres_seq(i, 1).name);
    end
end
