% compare_nonLTE_r49_versions
%
% Take nominal regr49 profiles and adjust solar zenith angles, then
%  compute TOA rads using different versions of SARTA built with different
%  non-LTE coefficients.
%  For use with CrIS HiRes.
%  Extended nonLTE solzen: 0    10    20    30    40    50    60    70    80    85    87    %       90    92    94  96    98   100   105   120
%
%
%
cd /home/chepplew/projects/sarta/prod_2018/cris_hr

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools            % subset_rtp
addpath /asl/matlib/aslutil             % int2bits
addpath /home/chepplew/myLib/matlib/rtptools    % replicate_rtp_headprof

% Setup temporary files
sTemp = mktemp();
fn_rtp1 = [sTemp '_' lower(csens) '_.rtp'];         % layers file
fn_rtp2 = [sTemp '_' lower(csens) '_sar.rtp'];      % calc file

% The nominal regr49 profiles:
rtp49 = '/home/chepplew/data/sarta/prod_2018/sarta_data/cris_hr/regr49_1013_400ppm.op.rtp';
[hd ha pd pa] = rtpread(rtp49);

% Setup for Sensor Specific. CrIS hiRes w/ 4 guard chans
load('/home/chepplew/projects/sarta/prod_2016/cris_hr/rtp_head_str_2235g4.mat');
hdr2.pfields = 3;
SARTAEXE = '/home/chepplew/gitLib/sarta/bin/crisg4_mar18_basic_optr_co2';

% For AIRS_L1C
hinfo = hdfinfo('/home/chepplew/projects/airs/L1C/airs_l1c_srf_tables_lls_new.hdf');
vchan = hdfread(hinfo.SDS(2));
ichan = hdfread(hinfo.SDS(1));
SARTAEXE = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_p2022jul22_v2';
SARTAEXE = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_p2022jul22_xnte';

% duplicate original profiles 8 times and add range of angles
satzen_angs = [0.0  8.8322  17.9223  32.8244  44.8285  53.4704  59.8336  65.0428];
solzen_angs = [22.5 22.5    22.5     22.5     22.5     22.5     22.5     22.5];
solzen_angs = [0.0  22.5    45.0     67.0     75.0     85.0     95.0     115.0];
solzen_angs = [0.0  10.0  20.0 30.0 40.0 50.0 60.0 70.0 80.0  85.0 87.0 ...
    90.0 92.0 94.0  98.0 100.0 105.0 120.0 150.0];
satzen_angs = [0.0  0.0 0.0 8.83 8.83 8.83 17.92 17.92 17.92 32.82 32.84 ... 
    44.82 44.82 53.47 53.47  59.83 59.83  65.04 65.04];
nangs = length(solzen_angs);

h = []; p = [];
for ii = 1 : 49
   [hy, py]  = replicate_rtp_headprof(hd,pd,ii,nangs);
   py.satzen = satzen_angs;
   py.solzen = solzen_angs;
   if ii == 1
     h = hy;
     p = py;
   else
     [h,p] = cat_rtp(h,p,hy,py);
   end
end

head = h;
head.nchan = length(vchan);
head.vchan = vchan;
head.ichan = ichan;

% Set up the SARTA run
rtpwrite(fn_rtp1, head, ha, p, pa);

unix([SARTAEXE ' fin=' fn_rtp1 ' fout=' fn_rtp2 ' > /home/chepplew/logs/sarta/sar_out.log']);

% Load results into memory
[hds has pds pas] = rtpread(fn_rtp2);

% Get angle subsets
idx.sol = [];
for i = 1:length(solzen_angs)
  idx.sol{i} = find(pds.solzen == solzen_angs(i));
end


% Calc BT and compute basic stats
btc = rad2bt(hds.vchan, pds.rcalc);
btc_mn = nanmean(btc,2);
btc_sd = nanstd(btc,0,2);

[sfreq, iss] = sort(hds.vchan);

%{
figure(1);clf;yyaxis left;plot(hds.vchan,btc_mn,'-');grid on;
  yyaxis right;plot(hds.vchan,btc_sd,'-');ylim([0 50]);xlim([640 2600])


%}
