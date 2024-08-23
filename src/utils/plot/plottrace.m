function fig = plottrace(pvt_seq)
    fig = figure('Position', [100, 100, 1024, 384]);
    subplot(121);
    pos_enu = [pvt_seq.PosENU];
    scatter(pos_enu(1, :), pos_enu(2, :), 10, 'filled'); grid on;
    xlabel("East (m)", "fontsize", 12)
    ylabel("North (m)", "fontsize", 12)
    hold on;
    
    % Drop positioning results with unacceptable error.
    pos_enu_f = pos_enu;
    pos_enu_f(:, abs(pos_enu(3, :))>20) = NaN;
    scatter(pos_enu_f(1, :), pos_enu_f(2, :), 10, 'filled');

    %% Time sequential analysis
    t = [pvt_seq.Time];
    v = [pvt_seq.Vel];
    rec_t_utc_sec = [3600, 60, 1]*[t.Hour; t.Minute; t.Second];
    dt = rec_t_utc_sec - rec_t_utc_sec(1);

    subplot(222);
    plot(dt, pos_enu(3, :),   'LineWidth', 1.0); grid on; hold on;
    plot(dt, pos_enu_f(3, :), 'LineWidth', 1.0); grid on
    xlabel("$t$(s)", "fontsize", 12, "Interpreter", "latex");
    ylabel("Altitude(m)", "fontsize", 12, "Interpreter", "latex");
    
    subplot(224);
    plot(dt, vecnorm(v), 'LineWidth', 1.0); grid on
    xlabel("$t$(s)", "fontsize", 12, "Interpreter", "latex");
    ylabel("Speed(m/s)", "fontSize", 12, "Interpreter", "latex");
end