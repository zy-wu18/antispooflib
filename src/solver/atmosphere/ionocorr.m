function [dIons, rhos_corr, drhos_corr, ps_corr, vs_corr, dts_corr, cnrs_corr, uobs_out] = ...
    ionocorr(ps, vs, dts, rhos, drhos, cnrs, alpha, beta, obs, uobs, upvt, opt, eph_dict)
switch opt
    case 'None'
        dIons = zeros(32, 1);
        rhos_corr = rhos;
        drhos_corr = drhos;
        ps_corr = ps;
        vs_corr = vs;
        dts_corr = dts;
        cnrs_corr = cnrs;
        uobs_out = uobs;
    case 'Klobuchar'
        
        upvt_LLA0 = upvt.PosLLA;
        az = [uobs.Az];
        el = [uobs.El];
        drhos_corr = drhos;
        [rhos_corr, dIons] = Klobuchar(alpha, beta, rhos, uobs, upvt_LLA0, az, el, eph_dict);
        ps_corr = ps;
        vs_corr = vs;
        dts_corr = dts;
        cnrs_corr = cnrs;
        uobs_out = uobs;
    otherwise
        [rhos_corr, drhos_corr, ps_corr, vs_corr, dts_corr, cnrs_corr, dIons, uobs_out] = DoubleFreq(ps, vs, dts, drhos, cnrs, obs, uobs, opt);
end