function obs_seq = readnvlog(fname)
    flines = readlines(fname);
    
    nobs_max = 86400;
    block_order = [1, nan, nan, nan, 3, nan]; % B1I first, then B3I.
    block_sync = false;
    
    obs_t = struct('Time', NaT, 'Sys', '?', 'PRN', NaN, 'Fc', NaN, ... 
        'Rho', NaN, 'ObsTime', NaN, 'Fd', NaN, 'AcPh', NaN, 'CNR', NaN);
    obs_seq = cell(1, nobs_max);
    nobs = 1;
    i = 1;

    while i < length(flines) && nobs <= nobs_max
        b = 1; % b = #obsevable_blocks
        while b <= length(block_order) && i < length(flines)
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
                
                if ~block_sync
                    if vals(field2idx("sig")) == block_order(b)
                        block_sync = true;
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
                    obs_seq{nobs} = insertobs(obs_seq{nobs}, obs); % append obs_seq
                end
            end
    
            if(block_sync)
                nobs = nobs + (b==length(block_order));
                b = b + 1; % next block
            end
        end
    end
    obs_seq = obs_seq(1:(nobs-1));
end

%% Parse NaVitech 'sig'
function [sys, fc] = sig2sysfc(sig)
    switch(sig)
        case 1
            sys = 'B';
            fc = 1561.098e6;
        case 3
            sys = 'B';
            fc = 1268.52e6;
        case 5
            sys = 'G';
            fc = 1575.42e6;
        case 6
            sys = 'B';
            fc = 1575.42e6;
        otherwise
            sys = '?';
            fc = nan;
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