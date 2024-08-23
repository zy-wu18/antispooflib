clearvars *;
addpath(genpath('utils\'));
addpath(genpath('reader\'));
addpath(genpath('solver\'));

try
    %% Load observation
    [obsfname, obsfpath] = uigetfile('*.obs', "Choose RINEX observable file...", "../data/");
    obs_seq = readrnx303([obsfpath, obsfname]);
    if(isempty(obs_seq))
        error("InvalidObservables: %s", [obsfpath, obsfname]);
    end
    
    %% Load broadcast ephemeris
    [ephfname, ephfpath] = uigetfile({ ...
            '*.2?n; *.2?g; *.rnx; *.nav', ...
            'RINEX304(*.2?n, *.2?g, *.rnx, *.nav)'}, ...
        "Choose RINEX304 emphemeris file...", "../data/");
    eph_dict = readrnx304([ephfpath, ephfname], obs_seq{1}(1).Time);

    %% Launch user PVT solver(1. solve P/V of satellites, 2. solve user P/V)
    pntcfg = struct( ... % PNT configurator
        'cnrMask', 20, 'elvMask', 10, ...
        'userLLA0', [30 120 0], ...
        'constellation', 'GRJEC', ...
        'pntSolver', @(rho, drho, ps, vs, cnrs)lse4pnt(rho, drho, ps, vs));
    ascfg = struct(... anti-spoofing configurator
        'spfDetector', [ ...
            struct('cfg', struct('NumAccum',30,'NumPoll',5,'Step',5,'Thresh',0.5), ...
                'detector', @(rho,drho,ps,vs,cnrs,keys,cfg)cnrcorr(cnrs,keys,cfg)), ...
            ...struct('cfg', struct('NumAccum',4,'NumPoll',5,'Pfa',1e-3), ...
            ...    'detector', @(rho,drho,ps,vs,cnrs,keys,cfg)draim2nd(drho,ps,vs,cfg)), ...
            struct('cfg', struct('NumAccum',1,'NumPoll',5,'Pfa',1e-3,'Sigma2',300), ...
                'detector', @(rho,drho,ps,vs,cnrs,keys,cfg)raim(rho,ps,cfg))], ...
        'spfRecognizer', []);
    [upvt_seq, uobs_seq, asres_seq] = launchpnt(obs_seq, eph_dict, pntcfg, ascfg);

    %% Visualization
    close all;
    plotobs(obs_seq, uobs_seq); % CNR, #SATV, #SATU plot
    plottrace(upvt_seq); % ENU trace
    plotspfdetect(asres_seq, [upvt_seq.Time]); % Spoof detection

catch WE
    if(strcmp(WE.identifier, 'MATLAB:fopen:NotScalarIntValued') || ...
       strcmp(WE.identifier, 'MATLAB:validators:mustBeNonzeroLengthText') || ...
       strcmp(WE.identifier, 'readrnx304:InvalidEph')|| ...
       strcmp(WE.identifier, 'MATLAB:notExistentField'))
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
