function obs_seq = readnvlog(fname)
    %% Initialization
    logger = Logger();
    logger.enStack("readnvlog: Loading from %s.", fname);
    flines = readlines(fname);
    
    nobs_max = 86400;
    blk_sig = [1, nan, nan, nan, 3, nan]; % blk_sig(i) = signal on the ith block
    blk_sync = false; % whether the observables of 1st block is found
    
    obs_t = struct('Time', NaT, 'Sys', '?', 'PRN', NaN, 'Fc', NaN, ... 
        'Rho', NaN, 'ObsTime', NaN, 'Fd', NaN, 'AcPh', NaN, 'CNR', NaN);
    obs_seq = cell(1, nobs_max);
    n = 1;
    i = 1;

    %% File -> DBCHH blocks (# = nobs) -> DBCHN lines (# = mobs)
    while i < length(flines) && n <= nobs_max
        logger.refreshBar(i, length(flines));
        b = 1; % reading the (b)th DBCH block of the (n)th frame
        while b <= length(blk_sig) && i < length(flines)
            fline = regexprep(flines{i}, '\s', '');
            i = i + 1;
            
            % $DBCHH: read channel header, start of a new block
            lbias_chh = strfind(fline, '$DBCHH');
            lformat = [];
            mobs = [];
            field2idx = dictionary();
        
            if(isempty(lbias_chh) && isempty(lformat)) % no history, no update
                continue;
            elseif(~isempty(lbias_chh)) % update lformat
                fline = fline((lbias_chh+7):end);
                lfields = strsplit(fline, {','});
                
                for j = 1:length(lfields)
                    lfield = lfields{j};
                    switch char(lfield)
                        case {'sig','ant','ch','sv','tk_sec','lk','currt','cut','iCorr','qCorr','eph','week'}
                            lformat = [lformat, '%d']; %#ok<*SAGROW,*AGROW> 
                        case {'aga','cn0','doppler','icv','agi','agq','cn0-i','sqmp','sqmd'}
                            lformat = [lformat, '%f'];
                        case {'transmittime','pseudo-range','carrierphase'}
                            lformat = [lformat, '%lf'];
                        otherwise
                            if(~isempty(str2double(lfield)) && ~isnan(str2double(lfield)) && isempty(mobs))
                                mobs = str2double(lfield); % first numeric field
                            end
                            lformat = [lformat, '%lf'];
                    end
                    field2idx(lfields{j}) = j;
                    if(j < length(lfields))
                        lformat = [lformat, ','];
                    end
                end
                assert(~isempty(mobs)); % the first numeric field in $DBCHH = #obs in this block
            end
    
            % $DBCHN: read channel information
            for j = 1:mobs
                fline = regexprep(flines{i}, '\s', '');
                i = i + 1;
                lbias_chn = strfind(fline, '$DBCHN');
                assert(~isempty(lbias_chn));
                assert(~isempty(lformat));
                fline = fline((lbias_chn+7):end);
                vals = sscanf(fline, lformat);
                
                if ~blk_sync % try block synchronization
                    if vals(field2idx("sig")) == blk_sig(1)
                        blk_sync = true;
                    else
                        continue;
                    end
                end
                
                obs = obs_t;
                [obs.Sys, obs.Fc] = sig2sysfc(vals(field2idx("sig")));
                obs.PRN = vals(field2idx("sv"));
                obs.Rho = vals(field2idx("pseudo-range"));
                obs.AcPh = vals(field2idx("carrierphase"));
                obs.Fd = vals(field2idx("doppler"));
                obs.CNR = vals(field2idx("cn0"));
                if(~isnan(obs.Rho) && obs.Rho>1 && ~isnan(obs.Sys))
                    obs_seq{n} = insertobs(obs_seq{n}, obs);
                end
            end
            
            if(blk_sync)
                if(b == length(blk_sig))
                    i = i + 1; % skip 1 line $DBANT
                    fline = regexprep(flines{i}, '\s', '');
                    lbias_chn = strfind(fline, '$GNZDA');
                    assert(~isempty(lbias_chn));
                    fline = fline((lbias_chn+7):end);
                    vals = sscanf(fline, '%02d%02d%02d.%02d,%d,%d,%d')';
                    rec_t_utc = vals([7,6,5,1,2,3])+[zeros(1,5), vals(4)/100];
                    i = i + 1;
                    fline = regexprep(flines{i}, '\s', '');
                    lbias_chn = strfind(fline, '$GNGGA');
                    assert(~isempty(lbias_chn));
                    fline = fline((lbias_chn+7):end);
                    vals = sscanf(fline, '%02d%02d%02d.%02d,%lf,%c,%lf,%c,%d,%d,%f,%f,%c,%f,%c');
                    if(rec_t_utc(1) >= 1980 && rec_t_utc(1) <= 2099)
                        rec_t_gps = Utc2Gps(rec_t_utc);
                        %rec_t_gps_wn = rec_t_gps(1);
                        rec_t_gps_tow  = rec_t_gps(2) - LeapSeconds(rec_t_utc); % Leap seconds after year 2017
                    else
                        rec_t_gps_tow = nan;
                    end
                    for j = 1:length(obs_seq{n})
                        obs_seq{n}(j).Time = rec_t_utc;
                        obs_seq{n}(j).ObsTime = rec_t_gps_tow;
                    end
                    n = n + 1; % next frame
                end
                b = b + 1; % next block            
            end
        end
    end
    obs_seq = obs_seq(1:(n-1));
    logger.resetBar;
    logger.writeLine("%d data blocks have been read successfully.", n-1);
    logger.writeLine("Recorded from %s to %s;", datetime(obs_seq{1}(1).Time), datetime(obs_seq{end}(1).Time));
    logger.writeLine("Maximum/Minimum #obs = %3d/%3d;", ...
        max([cellfun(@(x) length(x), obs_seq)]), ...
        min([cellfun(@(x) length(x), obs_seq)]));
    logger.deStack("readnvlog: %d observations loaded.\n", n-1);
end


%% Utils Parse NaVitech 'sig'
function [sys, fc] = sig2sysfc(sig)
    switch(sig)
        case 1 % B1I
            sys = 'B'; fc = 1561.098e6;
        case 3 % B3I
            sys = 'B'; fc = 1268.52e6;
        case 5 % GPS L1 C/A
            sys = 'G'; fc = 1575.42e6;
        case 6 % B1C
            sys = 'B'; fc = 1575.42e6;
        otherwise % undefined
            sys = '?'; fc = nan;
    end
end

%% Insert B3I after corresponding B1I if exists.
function obss = insertobs(obss, obs)
    if isempty(obss)
        obss = obs;
    else
        keys = arrayfun(@(a,b)sprintf("%c%02d",a,b), [obss.Sys], [obss.PRN]);
        query = sprintf("%c%02d", obs.Sys, obs.PRN);
        idx = find(arrayfun(@(a)strcmp(a, query), keys), 1);
        if(~isempty(idx))
            obss = [obss(1:idx(1)), obs, obss(idx(1)+1:end)];
        else
            obss = [obss, obs];
        end
    end
end