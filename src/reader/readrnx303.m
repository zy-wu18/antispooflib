function obs_seq = readrnx303(fname)
    logger = Logger();
    logger.enStack("readrnx303: Loading from %s.", fname);

    %% Output initialization
    ftext = fileread(fname); % Read file content as characters
    L = sum((ftext == '>')); % Number of blocks
    obs_t = struct('Time', NaT, 'Sys', '?', 'PRN', NaN, 'Fc', NaN, ... 
        'Rho', NaN, 'ObsTime', NaN, ...
        'Fd', NaN, 'AcPh', NaN, 'CNR', NaN);
    obs_seq = cell(1, L);
    logger.writeLine("%d data blocks detected.", L);
    
    %% Skip the header of .obs file
    idx2sname = struct('G', [], 'R', [], 'E', [], 'J', [], 'C', [], 'S', []);
    idx2field= struct('G', [], 'R', [], 'E', [], 'J', [], 'C', [], 'S', []);
    abbr2field = dictionary('C', 'Rho', 'L', 'AcPh', 'D', 'Fd', 'S', 'CNR');
    sys2nobs = dictionary('G', 0, 'R', 0, 'E', 0, 'J', 0, 'C', 0, 'S', 0);
    glonass_slot = zeros(1, 24);

    fid = fopen(fname);
    line = fgetl(fid);
    while ~contains(line, 'END OF HEADER')
        line = fgetl(fid);
        assert(ischar(line));
        if(contains(line, 'SYS / # / OBS TYPES'))
            if(line(1) ~= ' ') % only read the first 16 observations
                v = sscanf(line, '%c%d');
                sys = char(v(1));
                ssegs = strsplit(regexprep(line(7:60), '^\s+|\s+$', ''));
                assert(v(2) == length(ssegs));
                sys2nobs(sys) = v(2);
                for i = 1:v(2)
                    sseg = char(ssegs(i)); % e.g. 'S5Q'= strength of 5Q;
                    sname = [sys, sseg(2:3)]; % e.g. ['J', '5Q']->QZSS
                    idx2sname.(sys) = [idx2sname.(sys), string(sname)];
                    idx2field.(sys)= [idx2field.(sys), abbr2field(sseg(1))];
                end
            end
        elseif(contains(line, 'GLONASS SLOT / FRQ #'))
            line = regexprep(line(1:60), 'R', '');
            v = sscanf(line, '%d');
            glonass_slot(v(2:2:end)) = v(3:2:end);
        end
    end
    
    %% Read data, block by block
    idx = 0;
    while ~feof(fid)
        logger.refreshBar(idx, L);
        line = fgetl(fid);
        if line(1) == '>'
            idx = idx+1;
            line = line(2:end);
            rec_t_utc = [sscanf(line, "%d%d%d%d%d%f", 6)]';
            rec_t_gps = Utc2Gps(rec_t_utc);
            %rec_t_gps_wn = rec_t_gps(1);
            rec_t_gps_tow  = rec_t_gps(2) - LeapSeconds(rec_t_utc); % Leap seconds after year 2017
        elseif find(line(1) == ['G','R','E','J','C'], 1)
            obs = obs_t();
            obs.Time = rec_t_utc;
            obs.Sys = line(1);
            obs.PRN = str2double(line(2:3));
            for i = 1:sys2nobs(obs.Sys)
                sname = idx2sname.(obs.Sys)(i);
                obs.(idx2field.(obs.Sys)(i)) = str2double(line(4+((i*15-14):i*15)));
                obs.ObsTime = rec_t_gps_tow;
                obs.SigName = sname;
                
                if(i==sys2nobs(obs.Sys) || ~strcmp(sname, idx2sname.(obs.Sys)(i+1)))
                    obs.Fc = sname2fc(sname);
                    if(obs.Sys == 'R')
                        if(abs(obs.Fc-1602e6) < 1e-3)
                            obs.Fc = obs.Fc + glonass_slot(obs.PRN)*9e6/16;
                        elseif(abs(obs.Fc-1246e6) < 1e-3)
                            obs.Fc = obs.Fc + glonass_slot(obs.PRN)*7e6/16;
                        end
                    end
                    obs_seq{idx} = [obs_seq{idx}, obs]; % append obs_seq
                end
            end
        end
    end
    logger.resetBar;
    logger.writeLine("%d data blocks have been read successfully.", L);
    logger.writeLine("Recorded from %s to %s;", datetime(obs_seq{1}(1).Time), datetime(obs_seq{end}(1).Time));
    logger.writeLine("Maximum/Minimum #obs = %3d/%3d;", ...
        max([cellfun(@(x) length(x), obs_seq)]), ...
        min([cellfun(@(x) length(x), obs_seq)]));
    logger.deStack("readrnx303: %d observations loaded.\n", L);
    fclose(fid);
end

function fc_arr = sname2fc(sname)
% Frequency marks to center frequency of carriers.
    fc_arr = zeros(1, size(sname, 1));
    for i = 1:size(sname, 1)
        switch(sname)
            case {'G1A', 'G1B', 'G1C', 'G1P', 'G1W', 'E1A', 'E1B', 'E1C', 'E1X', 'E1Z', 'J1A', 'J1B', 'J1C', 'J1S'}
                fc = 1575.42e6;
            case {'G2C', 'G2P', 'G2W', 'G2S', 'G2L', 'G2X', 'J1L', 'J1X', 'J2S', 'J2L', 'J2X'}
                fc = 1227.6e6;
            case {'G5I', 'G5Q', 'G5X', 'E5I', 'E5Q', 'E5X', 'J5I', 'J5Q', 'J5X'}
                fc = 1176.45e6;
            case {'R1C', 'R1P'}
                fc = 1602e6;
            case {'R2C', 'R2P'}
                fc = 1246e6;
            case {'R3I', 'R3Q', 'R3X'}
                fc = 1202.025e6;
            case {'E7I', 'E7Q', 'E7X'}
                fc = 1207.140e6;
            case {'E8I', 'E8Q', 'E8X'}
                fc = 1191.795e6;
            case {'E6A', 'E6B', 'E6C', 'E6X', 'E6Z'}
                fc = 1278.75e6;
            case {'C2I', 'C2Q', 'C2X'}
                fc = 1561.098e6;
            case {'C7I', 'C7Q', 'C7X'}
                fc = 1207.14e6;
            case {'C6I', 'C6Q', 'C6X'}
                fc = 1268.52e6;
            otherwise 
                fc = NaN;
        end
        fc_arr(i) = fc;
    end
end
