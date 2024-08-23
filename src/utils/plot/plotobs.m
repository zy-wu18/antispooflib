function fig_obs = plotobs(obs_seq, uobs_seq)
% Plot observables
% args  :   1xL obs_t   obs_seq         full observable sequence
%           1xL obs_t   uobs_seq    [optional] observables used in PNT
% return:   figure      fig_obs
    if(nargin == 1)
        uobs_seq = obs_seq;
    end

    fig_obs = figure('Position', [110, 110, 896, 640]);
    L = size(obs_seq, 2);
    logger = Logger;
    logger.enStack("plotobs: Ploting observables.");
    
    %% Avalability view
    num_dict = dictionary();
    cnr = zeros(256, L)+NaN; % Assumed max #channel
    az = zeros(256, L)+NaN;
    el = zeros(256, L)+NaN;
    nsatv = zeros(5, L);
    nsatu = zeros(5, L);
    sys2num = dictionary('G',1,'R',2,'J',3,'E',4,'C',5);
    
    for i = 1:L
        obs = obs_seq{i};
        uobs = uobs_seq{i};
        M = length(obs);
        ju = 1;
        for j = 1:M
            o = obs(j);
            if(isempty(uobs))
                if(j > 1 && (o.PRN ~= obs(j-1).PRN || o.Sys ~= obs(j-1).Sys))
                    nsatv(sys2num(o.Sys), i) = nsatv(sys2num(o.Sys), i)+1;
                end
                continue;
            end
            ou = uobs(ju);
            
            % Skip unused obs, only +1 for the nsatv counter if not rebundant
            if(~isKey(sys2num, o.Sys) || ~strcmp(o.SigName, ou.SigName) || o.PRN ~= ou.PRN)
                if(j > 1 && (o.PRN ~= obs(j-1).PRN || o.Sys ~= obs(j-1).Sys))
                    nsatv(sys2num(o.Sys), i) = nsatv(sys2num(o.Sys), i)+1;
                end
                continue;
            end
            
            ju = min(length(uobs), ju + 1); % the order of observables are assumed to be kept
            nsatv(sys2num(o.Sys), i) = nsatv(sys2num(o.Sys), i)+1;
            nsatu(sys2num(o.Sys), i) = nsatu(sys2num(o.Sys), i)+1;
            
            key = sprintf("%c%02d", o.Sys, o.PRN);
            if(num_dict.numEntries == 0 || ~isKey(num_dict, key))
                % Insert dictionary
                num_dict(key) = num_dict.numEntries + 1;
            end
            cnr(num_dict(key), i) = ou.CNR;
            az(num_dict(key), i) = ou.Az;
            el(num_dict(key), i) = ou.El;
        end
    end
    nch = num_dict.numEntries;
    logger.writeLine("#detectedch=%d", nch);
    if(~isempty(uobs_seq{1}) && ~isempty(uobs_seq{end}))
        tstart = uobs_seq{1}(1).SatTime;
        tend = uobs_seq{end}(1).SatTime;
        logger.writeLine("tstart=%.1f[TOW], tend=%.1f[TOW]", tstart, tend);
        logger.writeLine("duration=%.1f[sec], L=%d", tend-tstart, L);
    else
        logger.writeLine("plotobs: uobs_seq start or end time error");
    end
    cnr = cnr(1:nch, :);
    az = az(1:nch, :);
    el = el(1:nch, :);
    
    %% CNR vision
    subplot(4,3,[1,2,4,5]);
        cnr_min = min(min(cnr(cnr~=0)));
        cnr_max = max(max(cnr(cnr~=0)));
        cmap = @(x)([1.0;0;0]+[-0.5;1;0]*max(0,(x-cnr_min)./(cnr_max-cnr_min),"includenan"));
        for ch = 1:nch
            logger.refreshBar(ch, nch);
            cnr_ch = cnr(ch, :);
            cnr_ch(cnr_ch==0) = NaN;
            scatter(1:L, repmat(ch, [1,L]), 5, cmap(cnr_ch)', 'filled', 'Marker', 'square');
            text(-L/10, ch, num_dict.keys{ch});
            hold on;
        end
        logger.resetBar;
        axis([0, L, 0, nch+1]); yticks([]);
    
    %% #SatUsed, #SatVis vision
    subplot(4,3,[7,8]);
        cssatv = cumsum(nsatv);
        plot(cssatv', 'LineWidth', 1); xlim([0, L]); ylim([0, max(cssatv(end, :))]);
        lgd = legend("#GPS", "+#GLO", "+#QZSS", "+#GAL", "+#BDS");
        set(lgd, "Location", "southeast");
        ylabel("#Visible sat"); xticks([]);
    subplot(4,3,[10,11]); 
        cssatu = cumsum(nsatu);
        plot(cssatu', 'LineWidth', 1); xlim([0, L]); ylim([0, max(cssatv(end, :))]);
        lgd = legend("#GPS", "+#GLO", "+#QZSS", "+#GAL", "+#BDS");
        set(lgd, "Location", "southeast");
        ylabel("#Used sat");
    
    %% Elevation and Azimuth angle -> skyplot vision
    subplot(2,3,3);
        marks = cellfun(@(s,p)sprintf("%c%02d",s,p), {uobs_seq{1}.Sys}, {uobs_seq{1}.PRN});
        az_d = az(~isnan(az(:, 1)), 1)/pi*180;
        el_d = el(~isnan(el(:, 1)), 1)/pi*180;
        grps = categorical({uobs_seq{1}.Sys}, {'G','R','J','S','E','C'});
        if(~isempty(az_d))
            skyplot(az_d, el_d, marks, 'GroupData', grps, 'MarkerSizeData', 15, 'LabelFontSize',8);
        end
    subplot(2,3,6);
        marks = cellfun(@(s,p)sprintf("%c%02d",s,p), {uobs_seq{end}.Sys}, {uobs_seq{end}.PRN});
        az_d = az(~isnan(az(:, end)), end)/pi*180;
        el_d = el(~isnan(el(:, end)), end)/pi*180;
        grps = categorical({uobs_seq{end}.Sys}, {'G','R','J','S','E','C'});
        if(~isempty(az_d))
            skyplot(az_d, el_d, marks, 'GroupData', grps, 'MarkerSizeData', 15, 'LabelFontSize',8);
        end
    logger.deStack("plotobs: finished plotting observables.");
end
