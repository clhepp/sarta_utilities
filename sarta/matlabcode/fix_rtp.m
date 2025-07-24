% fix the guard channels on the CrIS hiRes RTP file:

% from the following RTP file the gurad chanels are wrong and cause the sarta rtpopn.f
%  to fail.
%  max(hd.ichan) = 2233 (should be 2211 + 12 = 2223).
addpath /home/chepplew/myLib/matlib

run paths

dp = '/asl/s1/chepplew/projects/sarta/cris_hr/';
fn = 'rtp_d20150531_gran1_clear_hrg4.rtp';
fn = 'cris_hrg4_2015d151t1032228_clr.rtp';
[hd ha pd pa] = rtpread(strcat(dp,fn));

grdchns = find(hd.ichan > 2211);

% re-order these guard channels
%hd.ichan(grdchns)               % 2214:2217, 2222:2225, 2230:2233.
%hd.vchan(grdchns)
hd.ichan(grdchns(1:12)) = [2212:2223];
hd.ichan = sort(hd.ichan);
hd.vchan = sort(hd.vchan);

%{
figure(3);plot(hd.ichan,hd.vchan,'.');grid on; %xlim([2200 2250])
figure(3);plot(sort(hd.ichan), sort(hd.vchan),'.');grid on; %xlim([1 20])

%}

fo = 'rtp_d20150531_t0008263_clear_fix.rtp';
rtpwrite(strcat(dp,fo), hd, ha, pd, pa);

% OPTION 2
% re-arrange and make 4 guard channels per edge - adding new prof.robs1

dp = '/asl/rtp/rtp_cris_ccast_hires/clear/2015/151/';
fn = 'rtp_d20150531_t0008263_g4_clear.rtp';

[hd ha pd pa] = rtpread(strcat(dp,fn));
grdchns = find(hd.ichan > 2211);        % e.g. 1:4,718:725, 1591:1598, 2232:2235
[IX IY]  = sort(hd.ichan);
hd.ichan = sort(hd.ichan);
% hd.vchan are already in order of increasing wavenumber
hd.vchan = sort(hd.vchan);
% the robs1 appear to be aligned with frequency but some values are NAN


dout = '/asl/s1/chepplew/projects/sarta/cris_hr/';
[PATHSTR,NAME,EXT] = fileparts(fn);
fo   = [NAME '_fix.rtp'];

rtpwrite(strcat(dout,fo), hd, ha, pd, pa);

