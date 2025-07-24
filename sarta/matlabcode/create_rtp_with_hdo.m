% create_rtp_with_hdo
%
% Takes an existing rtp clear subset file, adds HDO depletion to prof.uded(20) 
% One value of depletion per profile
% For use with modified kCARTa and SARTA
%
%


addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath /asl/matlib/plotutils                   % aslprint
addpath /asl/packages/airs_decon/source         % hamm_app
addpath /home/sbuczko1/git/rtp_prod2/util       % rtp_sub_prof

% Choose sensors
csens = 'AIRS_L1C';

% KLAYERS
klayersexe = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';


% Get spectral grid to update head structure and the two sarta execs
switch csens
  case 'AIRS_L1C'
    fn_rtp1 = '/asl/rtp/airs/airs_l1c_v674/clear/2021/ecmwf_airicrad_day270_clear.rtp';
    vchan = hdfread('/home/chepplew/projects/airs/L1C/airs_l1c_srf_tables_lls_new.hdf','freq');
    ichan = hdfread('/home/chepplew/projects/airs/L1C/airs_l1c_srf_tables_lls_new.hdf','chanid');
  case 'CRIS_FSR'
    fn_rtp1 = '/asl/rtp/cris/j01_ccast_hires/clear/2019/cris2_era_csarta_clear_d20190621.rtp';

  case 'CRIS_NSR'

  case 'CHIRP'

  case 'IASI'

end

% Get the original clear subset RTP file:
[head hatt prof patt] = rtpread(fn_rtp1);

% reduce number of profiles down
iidx = [1:length(prof.rlat)];
jidx = sort(datasample(iidx,5000));

sub_prof = rtp_sub_prof(prof,jidx);



