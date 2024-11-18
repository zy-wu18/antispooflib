function [dions, rhos_corr, uobs] = ...
    ionocorr(rhos, obs, uobs, upvt, opt, eph_dict)
% Calculate ionospheric delay
% args  :   3xM double  ps      [m], satellite ECEF position [x, y, z]
%           3XM double  vs      [m/s], satellite ECEF velocity [x, y, z]
%           1xM double  dts     [s], satellite clock fix
%           1xM+ struct obs     satellites dual frequency observation, with
%                               2nd frequency and filtered ones included
%           1xM struct  uobs    satellites observation used in PVT calculation
%                               for IF, uobs only contains L1 observations.
%           1x1 struct  upvt    previous receiver PVT observation
%           strig       opt     Ionospheric correction model selection,
%                               'N' for None;'K' for Klobuchar model;
%                               'IF' for dual frequency model;
%           dictionary eph_dict ephemeris data struct, including ionospheric
%                               parameters.
% return:   1xM double  dions       [s], ionospheric delay
%           1xM double  rhos_corr   [m], satellite pseodoranges
%           1xM struct  uobs        satellites observation used in PVT calculation
    
    c = 2.99792458e8;
    M = length(uobs);
    dions = zeros(1, M);
    rhos_corr = zeros(1, M);
    
    switch upper(opt)
        case {'K', 'KLOBUCHAR'}
            [rhos_corr, dions] = Klobuchar(uobs, upvt.PosLLA, eph_dict);
        case {'IF', 'IFLC', 'IONOFREE'}
            [rhos_corr, dions, uobs] = dualfreq(obs, uobs, eph_dict);
        case 'GIM'
            ionpath = 'E:/Seafile/GNSS/gnss.workspace/ionstatistics/ionstatistics/full_ion_stat_2023-2024_1-302.mat';
            [rhos_corr, dions] = gim2dions(ionpath, uobs, upvt);
        case {'N', 'NULL'}
            for k = 1:M
                key = sprintf("%c%02d", uobs(k).Sys, uobs(k).PRN);
                if(eph_dict.isKey(key) && ~isempty(eph_dict(key).TGD))
                    TGD = eph_dict(key).TGD;
                    rhos_corr(k) = rhos(k) - c*TGD(1);
                else
                    rhos_corr(k) = rhos(k);
                end
            end
        otherwise
            for k = 1:M
                key = sprintf("%c%02d", uobs(k).Sys, uobs(k).PRN);
                if(eph_dict.isKey(key) && ~isempty(eph_dict(key).TGD))
                    TGD = eph_dict(key).TGD;
                    rhos_corr(k) = rhos(k) - c*TGD(1);
                else
                    rhos_corr(k) = rhos(k);
                end
            end
    end
end