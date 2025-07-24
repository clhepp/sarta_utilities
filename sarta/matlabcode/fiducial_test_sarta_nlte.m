% Day & Night Test of nonLTE using SARTA 
% for Sergio kCARTA comparison
% Jan 2020

addpath /asl/matlib/h4tools                        % rtpread
addpath /home/chepplew/projects/sarta/matlabcode
addpath /asl/matlib/rtptools/                      % rtpwrite_12

cd /home/chepplew/projects/sarta/prod_2019/

SARTAEXE = '/home/chepplew/gitLib/sarta/bin/iasi_jun19_co2nte400';

% Origin atmospheric test profiles: 
rtp.fin = ['/home/chepplew/data/sarta/prod_2019/generic/' ...
           'r49_1013_400p_emis1_7angs_unpert.rtp'];

[hd ha pd pa] = rtpread(rtp.fin);

    x = load('/home/chepplew/gitLib/ftc_dev/chanLists/iasi_8461_list');
    hd.vchan = single(x(:,2));
    hd.ichan = int32(x(:,1));
    hd.nchan = length(x);
    fnrtpx = ['outd/r49_1013m_400p_e1_7ang_unpert_8461.rtp'];
    tmpx = mktemp();
    outfiles = rtpwrite_12(tmpx,hd,ha,pd,pa);


