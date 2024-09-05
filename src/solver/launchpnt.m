function [upvt_seq, uobs_seq, asres_seq] = launchpnt(obs_seq, eph_dict, pntcfg, ascfg)
% launch PNT core with observables ranged in sequence and eph dictionary
% args  :   1xL cell    obs_seq     observables ranged in sequence, where
%                                   each element represents an obs array
%           dictionary  eph_dict    string->struct(neph_t, or geph_t) with
%                                   key=sprintf("%c%02d", sys, prn)
%           pntcfg_t    pntcfg      PNT configuration, includes
%                                   elvMask([deg], default=0),
%                                   constellation(char, default='GRJEC'), 
%                                   cnrMask([dBHz], default=0), 
%                                   pntSolver(handler@(rho, drho, ps, vs)),
%           ascfg_t     ascfg       anti-spoofing configuration, include
%                                   spoofDetector(sd_t, detector+cfg).
% return:   1xL pvt_t   upvt_seq    user's PVT solution result sequence
%           1xL cell    uobs_seq    used observables sequence, each element
%                                   would be 1xM obs_t
%           1xL asres_t asres_seq   anti-spoofing result sequence, includes
%                                   detection & recognition<TO DO>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(nargin <= 2)
        pntcfg = struct();
        pntcfg.elvMask = 0;
        pntcfg.constellation = 'GRJEC';
        pntcfg.cnrMask = 0;
        pntcfg.pntSolver = @(rho, drho, ps, vs, cnrs)lse4pnt(rho, drho, ps, vs);
        pntcfg.userLLA0 = [nan, nan, nan];
    elseif(nargin <= 3)
        ascfg.spfDetector = [];
        ascfg.spfRecognizer = [];
    end

    %% Initialize necessary parameters  
    L = size(obs_seq, 2);
    NaN3 = zeros(1, 3) + NaN;
    pvt_t = struct("Pos", NaN3, "Vel", NaN3, "Time", NaT, "Drift", NaN, ...
                   "PosENU", NaN3, "PosLLA", NaN3);
    asres_t = struct('name', [], 'alarm', [], 'T', [], 'gamma', []);
    
    logger = Logger();
    logger.enStack("launchpnt: launching PNT solution core, L=%d", L);
    logger.writeLine("pntcfg.constellation=%s", pntcfg.constellation);
    logger.writeLine("pntcfg.cnrMask=%.1f[dBHz]", pntcfg.cnrMask);
    logger.writeLine("pntcfg.elvMask=%.1f[deg]", pntcfg.elvMask);
    logger.writeLine("pntcfg.userLLA0=[%.6f,%.6f,%.1f]", pntcfg.userLLA0);
    logger.enStack("pntcfg.pntSolver=");
    logger.writeLine("%s",func2str(pntcfg.pntSolver));
    logger.deStack();
    if(~isempty(ascfg.spfDetector))
        logger.enStack("ascfg.spfDetector=(#%d)array", length(ascfg.spfDetector));
        cellfun(@(sd)logger.writeLine("%s",func2str(sd)), {ascfg.spfDetector(:).detector});
        logger.deStack();
    else
        logger.writeLine("ascfg.spfDetector=[]");
    end
    if(~isempty(ascfg.spfRecognizer))
        logger.enStack("ascfg.spfRecognizer=");
        cellfun(@(sd)logger.writeLine("%s",func2str(sd)), {ascfg.spfDetector(:).recognizer});
        logger.deStack('\n');
    else
        logger.writeLine("ascfg.spfRecognizer=[]\n");
    end
    
    
    %% Calculate position and velocity of satellites, then PVT solution 
    logger.enStack("pntSolver: User PNT calculator launched...");
    uobs_seq = cell([1, L]);
    upvt_seq = repmat(pvt_t, [1, L]);
    asres_seq = repmat(asres_t, [length(ascfg.spfDetector), L]);
    
    for i = 1:L
        logger.refreshBar(i, L);
        
        % Satellite position and observables fix
        [ps, vs, dts, rhos, drhos, cnrs, uobs_seq{i}] = ...
            ephposfix(obs_seq{i}, eph_dict, pntcfg, lla2ecef(pntcfg.userLLA0));
        % PNT solver
        [pu, vu, dtu, ddtu] = ...
            pntcfg.pntSolver(rhos', drhos', ps', vs', dts', cnrs', [uobs_seq{i}.Sys]);
        upvt_seq(i).Pos = pu';
        upvt_seq(i).Vel = vu';
        if(isempty(uobs_seq{i}))
            upvt_seq(i).Time = NaT;
        else
            upvt_seq(i).Time = datetime(uobs_seq{i}(1).Time) + dtu(1);
        end
        upvt_seq(i).Drift = ddtu;
        pntcfg.userLLA0 = ecef2lla(upvt_seq(max(i-1, 1)).Pos');
        
        % Spoofing detection result
        if(~isempty(ascfg.spfDetector))
            keys = arrayfun(@(c,p)sprintf("%c%02d",c,p), [uobs_seq{i}.Sys], [uobs_seq{i}.PRN]);
            sdres = spoofdetect(rhos',drhos',ps',vs',dts',cnrs',keys',ascfg.spfDetector);
            asres_seq(:, i) = sdres;
        end
    end
    logger.resetBar();
    logger.writeLine("User status obtained in ECEF coordinate.");

    %% ECEF to ENU
    lla0 = ecef2lla(upvt_seq(end).Pos');
    user_pos = [upvt_seq.Pos];
    [E, N, U] = ecef2enu(user_pos(1, :), user_pos(2, :), user_pos(3, :), ...
        lla0(1), lla0(2), lla0(3), wgs84Ellipsoid('meter'));
    for i = 1:L
        upvt_seq(i).PosLLA = ecef2lla(upvt_seq(i).Pos');
        upvt_seq(i).PosENU = [E(i); N(i); U(i)];
    end
    logger.writeLine("User status obtained in ENU coordinate.");
    logger.deStack("pntSolver: User PVT calculation finished.");

    logger.deStack("launchpnt: user PNT finished.\n");
end

