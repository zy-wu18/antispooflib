obs = obs_seq{1};
m = length(obs);
for i = 1:m
    sys = obs(i).Sys;
    prn = obs(i).PRN;
    rho = obs(i).Rho;
    acph = obs(i).AcPh;
    fd = obs(i).Fd;
    cnr = obs(i).CNR;
    if(strcmp(obs(i).SigName, 'C2I'))
        fprintf("%c%02d%14.3f%17.4f%15.3f%16.3f  ", sys, prn, rho, acph, fd, cnr);
        fprintf("\n");
    end
end