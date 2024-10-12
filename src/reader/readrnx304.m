function eph_dict = readrnx304(fname, t)
% RINEX304_Reader: yield ephemeris in <fname> nearest to datetime <t>
% args  :   string      fname   full path the RINEX 3.04 navigation data file
%           datetime    t       query UTC time
% return:   dictionary  eph     1x3 char -> neph_t or geph_t
% notes :   readrnx304 trys to gather all the ephemeris in <fname> whose
%           UTC time is close enough with <t>.
%           The closest epheremis of each satellite(e.g. 'G01') can be 
%           found in eph_dict('G01') if at least one ephemeris of G01 with 
%           a close-enough UTC time(<2 hours) was found in <fname>.

    %% Initialization
    logger = Logger();
    logger.enStack("readrnx304: Loading from %s.", fname);
    [~, ~, fext] = fileparts(fname);
    fid = fopen(fname);
    line = fgetl(fid);
    if(strcmp(fext, '.rnx') || contains(line, 'GNSS'))
        nav_fmt_l1 = "%c%d%d%d%d%d%d%d%lf%lf%lf"; % format of line 1
        nav_fmt_ln = "%lf%lf%lf%lf"; % format of line n > 1
        sys_default = '';
    elseif(fext(end) == 'n' || contains(line, 'GPS')) % GPS format
        nav_fmt_l1 = "%d%d%d%d%d%d%lf%lf%lf%lf"; % format of line 1
        nav_fmt_ln = "%lf%lf%lf%lf"; % format of line n > 1
        sys_default = 'G';
    elseif(fext(end) == 'g') % GLONASS format
        nav_fmt_l1 = "%d%d%d%d%d%d%lf%lf%lf%lf"; % format of line 1
        nav_fmt_ln = "%lf%lf%lf%lf"; % format of line n > 1
        sys_default = 'R';
    end
    
    %% Skip header
    line = fgetl(fid);
    l = 1;
    while ~contains(line, 'END OF HEADER')
        line = fgetl(fid);
        l = l+1;
    end
    
    %% Output initialization
    eph_dict = dictionary();
    
    % GPS, Galileo, BDS, QZNSS ephemeris
    eph_n_t = struct('sys',[],  'PRN',[],       'Toc_cal',[], ...
        'Toc_WN', [],           'Toc_TOW', [],  'Toes', [], ...
        'af0',[],   'af1',[],   'af2',[], ...
        'IODE',[],  'Crs',[],   'Delta_n',[],   'M_0',[],...
        'Cuc',[],   'e',[],     'Cus',[],       'sqrt_a',[], ...
        'Toe',[],   'Cic',[],   'Omega_0',[],   'Cis',[], ...
        'i_0',[],   'Crc',[],   'omega', [],    'Omega_dot',[], ...
        'I_dot',[], 'WN', [],   'Health', [],   'TGD',[]);
    
    % GLONASS, SBAS ephemeris
    eph_g_t = struct('sys', [], 'PRN',[],'Toc_cal',[], ...
        'Toc_WN',[], 'Toc_TOW', [], 'toc', [], ...
        'TauN',[], 'GammaN',[], 'tod',[], 'toe',[], ...
        'X', [], 'Xvel', [], 'Xacc', [], 'Health', [], ...
        'Y', [], 'Yvel', [], 'Yacc', [], 'k', [], ...
        'Z', [], 'Zvel', [], 'Zacc', [], 'age', []);
    
    %% Search for the nearest available ephemeris
    duration_tol = duration(47, 59, 59);
    delta_t = struct('G', repmat(duration_tol, [1, 32]), ...
        'E', repmat(duration_tol, [1, 36]), ...
        'C', repmat(duration_tol, [1, 64]), ...
        'J', repmat(duration_tol, [1,  7]), ...
        'R', repmat(duration_tol, [1, 26]), ...
        'S', repmat(duration_tol, [1, 60]));
    t = datetime(int32(t));
    num_sys = zeros(1, 6);
    logger.enStack("Searching closest ephemeris available, duration_tol=%s", duration_tol);
    logger.setTable(["%3d", "%3d", "%3d", "%3d", "%3d", "%3d"], ...
        ["G(GPS)", "C(BDS)", "J(QZSS)", "S(SBAS)", "E(Galileo)", "R(GLONASS)"]);

    while ~feof(fid)
        line = replace(fgetl(fid), "D", "e");
        [v, vn] = sscanf(line, nav_fmt_l1);
        if(any(char(v(1)) == ['G' 'E' 'C' 'J']) || (strcmp(sys_default,'G') && line(2)~=' '))
            eph_tmp = eph_n_t();
        elseif(any(char(v(1)) == ['R' 'S']) || (strcmp(sys_default,'R') && line(2)~=' '))
            eph_tmp = eph_g_t();
        else
            continue;
        end
        assert(vn + length(sys_default) == 11);
        if(isempty(sys_default))
            eph_tmp.sys = char(v(1));
        else
            eph_tmp.sys = sys_default;
        end
        eph_tmp.PRN = v(2-length(sys_default));
        v(3-length(sys_default)) = (v(3-length(sys_default)) < 100)*2000 + v(3-length(sys_default));
        eph_tmp.Toc_cal = datetime(v((3:8)-length(sys_default))');
        sys_time_bias = 0;
        if(eph_tmp.sys == 'C')
            sys_time_bias = 14;
        end
        eph_tmp.Toc_cal = eph_tmp.Toc_cal + duration(0, 0, sys_time_bias);
        gnss_time = Utc2Gps(v((3:8)-length(sys_default))'+[0 0 0 0 0 sys_time_bias]);
        eph_tmp.Toc_WN = gnss_time(1);
        eph_tmp.Toc_TOW= gnss_time(2);

        if(abs(duration(eph_tmp.Toc_cal - t)) < delta_t.(char(eph_tmp.sys))(eph_tmp.PRN))
            delta_t.(char(eph_tmp.sys))(eph_tmp.PRN) = abs(duration(eph_tmp.Toc_cal - t));
        else
            continue;
        end
        
        if(any(char(eph_tmp.sys) == ['G' 'E' 'C' 'J']))
            eph_tmp.af0 = v(9-length(sys_default));
            eph_tmp.af1 = v(10-length(sys_default));
            eph_tmp.af2 = v(11-length(sys_default));
            line = replace(fgetl(fid), "D", "e");
            [v, vn] = sscanf(line, nav_fmt_ln); assert(vn == 4);
            eph_tmp.IODE = v(1);
            eph_tmp.Crs = v(2);
            eph_tmp.Delta_n = v(3);
            eph_tmp.M_0 = v(4);
            line = replace(fgetl(fid), "D", "e");
            [v, vn] = sscanf(line, nav_fmt_ln); assert(vn == 4);
            eph_tmp.Cuc = v(1);
            eph_tmp.e = v(2);
            eph_tmp.Cus = v(3);
            eph_tmp.sqrt_a = v(4);
            line = replace(fgetl(fid), "D", "e");
            [v, vn] = sscanf(line, nav_fmt_ln); assert(vn == 4);
            eph_tmp.Toe = v(1) + sys_time_bias;
            eph_tmp.Toes = v(1);
            eph_tmp.Cic = v(2);
            eph_tmp.Omega_0 = v(3);
            eph_tmp.Cis = v(4);
            line = replace(fgetl(fid), "D", "e");
            [v, vn] = sscanf(line, nav_fmt_ln); assert(vn == 4);
            eph_tmp.i_0 = v(1);
            eph_tmp.Crc = v(2);
            eph_tmp.omega = v(3);
            eph_tmp.Omega_dot = v(4);
            line = replace(fgetl(fid), "D", "e");
            [v, vn] = sscanf(line, nav_fmt_ln); assert(vn >= 3);
            eph_tmp.I_dot = v(1);
            eph_tmp.WN = v(3);
            line = replace(fgetl(fid), "D", "e");
            [v, vn] = sscanf(line, nav_fmt_ln); assert(vn >= 3);
            if(eph_tmp.sys == 'J')
                eph_tmp.Health = bitand(v(2), 254);
            else
                eph_tmp.Health = v(2);
            end
            if(any(eph_tmp.sys == 'EC'))
                eph_tmp.TGD = [v(3), v(4)];
            else
                eph_tmp.TGD = v(3);
            end
            fgetl(fid);
        elseif(any(char(eph_tmp.sys) == ['R' 'S']))
            eph_tmp.TauN = -1.0*v(9-length(sys_default)); % ! -TauN in [sec]
            eph_tmp.GammaN = v(10-length(sys_default));
            eph_tmp.tod = v(11-length(sys_default));
            eph_tmp.toc = round15minutes(eph_tmp.Toc_cal);
            toc = eph_tmp.toc;
            eph_tmp.toe = Utc2Gps([toc.Year, toc.Month, toc.Day, toc.Hour, toc.Minute, toc.Second]);
            eph_tmp.toe = eph_tmp.toe(2);
            line = replace(fgetl(fid), "D", "e");
            [v, vn] = sscanf(line, nav_fmt_ln); assert(vn == 4);
            eph_tmp.X = v(1);
            eph_tmp.Xvel = v(2);
            eph_tmp.Xacc = v(3);
            eph_tmp.Health = v(4);
            line = replace(fgetl(fid), "D", "e");
            [v, vn] = sscanf(line, nav_fmt_ln); assert(vn == 4);
            eph_tmp.Y = v(1);
            eph_tmp.Yvel = v(2);
            eph_tmp.Yacc = v(3);
            eph_tmp.k = v(4);
            line = replace(fgetl(fid), "D", "e");
            [v, vn] = sscanf(line, nav_fmt_ln); assert(vn == 4);
            eph_tmp.Z = v(1);
            eph_tmp.Zvel = v(2);
            eph_tmp.Zacc = v(3);
            eph_tmp.age = v(4);
        end
        key = sprintf("%c%02d", eph_tmp.sys, eph_tmp.PRN);
        eph_dict(key) = eph_tmp;
        num_sys_last = num_sys;
        num_sys = arrayfun(@(sys)sum(startsWith(eph_dict.keys, sys)), 'GCJSER')';
        if(mod(sum(num_sys), 10) == 0 && any(num_sys_last ~= num_sys))
            logger.refreshTable(num_sys);
        end
    end
    fclose(fid);
    logger.deStack();
    assert(eph_dict.numEntries > 0, 'readrnx304:InvalidEph', ...
        sprintf("Ephemeris invalid:\n  nav time = %s;\n  obs time = %s.\n", eph_tmp.Toc_cal, t));

    logger.enStack("readrnx304 report:");
    for sys = ['G', 'C', 'J', 'S', 'E', 'R']
        num_sys = sum(startsWith(eph_dict.keys, sys));
        if(num_sys > 0)
            logger.writeLine("%c: %d ephemeris valid.", sys, num_sys);
        end
    end
    logger.deStack("Number of valid emphemeris: %d.", eph_dict.numEntries);
    logger.deStack("readrnx304 finished.\n");
end