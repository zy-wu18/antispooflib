clearvars *;
addpath(genpath('utils\'));
addpath(genpath('reader\'));
addpath(genpath('solver\'));

try
    %% Load observation
    [obsfname, obsfpath] = uigetfile({'*.obs'; '*.log'}, ...
        "Choose observable file or Navi-log file...", "../data/");
    [~, ~, obsfext] = fileparts(obsfname); 
    switch(obsfext)
        case '.obs'
            obs_seq = readobs([obsfpath, obsfname]);
        case '.log'
            obs_seq = readnvlog([obsfpath, obsfname]);
        otherwise
            error('Read observation error. Unsupported extention (%s).', obsfext);
    end
    if(isempty(obs_seq))
        error("InvalidObservables: %s", [obsfpath, obsfname]);
    end
    
    %% Load broadcast ephemeris
    [ephfname, ephfpath] = uigetfile({ ...
            '*.2?n; *.2?g; *.rnx; *.nav', ...
            'RINEX304(*.2?n, *.2?g, *.rnx, *.nav)'}, ...
        "Choose RINEX304 emphemeris file...", "../data/");
    eph_dict = readnav([ephfpath, ephfname], obs_seq{1}(1).Time);

    %% Launch user PVT solver(1. solve P/V of satellites, 2. solve user P/V)
    obs_rate = 10; % [Hz], observation rate
    pntcfg = struct( ... % PNT configurator
        'cnrMask', 0, 'elvMask', 0, ...
        'userLLA0', [40 116 0], ...
        'ionoOpt', 'IF', ... % N: Null; K: Klobuchar; IF: Iono-Free
        'constellation', 'GJEC', ...
        'pntSolver', @(rho,drho,ps,vs,dts,cnrs,sys)lse4pnt(rho,drho,ps,vs,dts,sys));
    ascfg = struct(... anti-spoofing configurator
        'spfDetector', [ ...
            struct('cfg', struct('ObsRate',obs_rate,'TAccum',3.0,'NumPoll',5,'TStep',1.0,'Thresh',0.5), ...
                'detector', @(rho,drho,ps,vs,dts,cnrs,keys,cfg)cnrcorr(cnrs,keys,cfg)), ...
            struct('cfg', struct('ObsRate',obs_rate,'TAccum',1.0,'NumPoll',5,'Pfa',1e-3), ...
                'detector', @(rho,drho,ps,vs,dts,cnrs,keys,cfg)draim2nd(drho,ps,vs,cfg)), ...
            struct('cfg', struct('ObsRate',obs_rate,'TAccum',1.0,'NumPoll',5,'Af2Thresh',5e-8), ...
                'detector', @(rho,drho,ps,vs,dts,cnrs,keys,cfg)cdm(drho,ps,vs,cfg)), ...
            struct('cfg', struct('ObsRate',obs_rate,'TAccum',0.1,'NumPoll',5,'Pfa',1e-3,'Sigma2',300), ...
                'detector', @(rho,drho,ps,vs,dts,cnrs,keys,cfg)raim(rho,ps,dts,keys,cfg))], ...
        'spfRecognizer', []);
    [upvt_seq, uobs_seq, asres_seq] = launchpnt(obs_seq, eph_dict, pntcfg, ascfg);

    %% Visualization
    close all;
    plotobs(obs_seq, uobs_seq); % CNR, #SATV, #SATU plot
    plottrace(upvt_seq); % ENU trace
    plotspfdetect(asres_seq, obs_rate); % Spoof detection

catch WE
    if(strcmp(WE.identifier, 'MATLAB:fopen:NotScalarIntValued') || ...
       strcmp(WE.identifier, 'MATLAB:validators:mustBeNonzeroLengthText') || ...
       strcmp(WE.identifier, 'readrnx304:InvalidEph')|| ...
       strcmp(WE.identifier, 'MATLAB:notExistentField')|| ...
       strcmp(WE.identifier, 'MATLAB:fileparts:MustBeChar'))
        fprintf("Known error encountered!\n  %s: %s.\n", WE.identifier, WE.message);
        arrayfun(@(x) fprintf("  %s, Line %d.\n", x.name, x.line), WE.stack); 
        fprintf("\n");
        fclose all;
        return;
    else
        fclose all;
        rethrow(WE);
    end
end
