function obs_seq = readnvlog(fname)
    %% Initialization
    logger = Logger();
    logger.enStack("readnvlog: Loading from %s.", fname);
    flines = readlines(fname);
    
    nobs_max = 86400;
    blk_sig = [1, nan, nan, nan, 3, nan]; % blk_sig(i) = signal on the ith block
    blk_sync = false; % whether the observables of 1st block is found
    
    obs_t = struct('Time', NaT, 'Sys', '?', 'PRN', NaN, 'SigName', [], ... 
        'ObsTime', NaN, 'Fc', NaN, 'Rho', NaN, 'Fd', NaN, 'AcPh', NaN, 'CNR', NaN);
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
                [obs.Sys, obs.Fc, obs.SigName] = sig2sysfc(vals(field2idx("sig")));
                obs.PRN = vals(field2idx("sv"));
                obs.Rho = vals(field2idx("pseudo-range"));
                obs.AcPh = vals(field2idx("carrierphase"));
                obs.Fd = vals(field2idx("doppler"));
                obs.CNR = vals(field2idx("cn0"));
                obs.ObsTime = vals(field2idx("transmittime")) + obs.Rho/299792458;
                if(obs.Sys == 'C')
                    obs.ObsTime = mod(obs.ObsTime + 14.0, 7*86400);
                end
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
                    yrs = vals(7);
                    mon = vals(6);
                    day = vals(5);
                    i = i + 1;
                    fline = regexprep(flines{i}, '\s', '');
                    lbias_chn = strfind(fline, '$GNGGA');
                    assert(~isempty(lbias_chn));
                    fline = fline((lbias_chn+7):end);
                    vals = sscanf(fline, '%02d%02d%02d.%02d,%lf,%c,%lf,%c,%d,%d,%f,%f,%c,%f,%c');
                    if(yrs >= 1980 && yrs <= 2099)
                        for j = 1:length(obs_seq{n})
                            dow = floor(obs.ObsTime/86400);
                            hrs = floor((obs.ObsTime - dow*86400)/3600);
                            mnt = floor((obs.ObsTime - (dow*86400+hrs*3600))/60);
                            sec = obs.ObsTime - ((dow*24+hrs)*60+mnt)*60;
                            obs_seq{n}(j).Time = [yrs, mon, day, hrs, mnt, sec];
                        end
                        n = n + 1; % next frame
                    else
                        obs_seq{n} = [];
                    end
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
function [sys, fc, signame] = sig2sysfc(sig)
    switch(sig)
        case 1 % B1I
            sys = 'C'; fc = 1561.098e6; signame = 'C2I';
        case 3 % B3I
            sys = 'C'; fc = 1268.52e6; signame = 'C6I';
        case 5 % GPS L1 C/A
            sys = 'G'; fc = 1575.42e6; signame = 'G1C';
        case 6 % B1C
            sys = 'C'; fc = 1561.098e6; signame = 'C2X';
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