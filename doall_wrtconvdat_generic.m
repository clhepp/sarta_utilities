function [ok]=doall_wrtconvdat_generic(opts, comment, pnums, lfcow2fow);

% function [ok]=doall_wrtconvdat(comment, pnums, rtpfile, matfile, outfile,lfcow2fow);
%
% Merge the matlab convolved l-to-s trans data from kcarta with the fortran profile
% data and write out a fortran binary data.
% One output file per profile. Used by fitftc.
%
% Input: opts [ structure]
%       csens    = {string} ('cris', 'airs', 'iasi' 'chirp') the sensor to use.
%       regset = [string] ('r49','saf704') which regression profile set to use.
%       myset = {string} set1 ... set5 (set 5 is used for 5, 6 and 7)
%       nscang: number of scan angles used by kCARTA in the OD L2S files.
%
%    comment  = {string (max 35 char)}, comment string which will be
%       Example: ['6 Aug 99, SRF 149 Aug99, fow, p2019']
%    pnums = [1 x nprof] integer number of profiles {eg [1:48]}
%    lfcow2fow = OPTIONAL [1 x 1] convert fcowB3 input to fowB3
%       output? (0=no=default, 1=yes). When "yes", comment and outfile
%       should be for fowB3.
%
% Other local params set in script:
%    rtpfile = {string} name of RTP file with regression data profile
%    matfile = {string} name prefix of matlab convolved trans files.
%       The file names are assumed to be:
%       <whatever><pnum>'.mat'
%    outfile = string, name prefix of fortran output file to create.
%       The output file names will be:
%       <outfile><pnum>'.dat'
%
% Output:
%    ok = integer, 0 if the function detects a problem, otherwise 1
%
% Child procs: wrt_convdat.m
%
% Note on ordering layers:
%   The convolved kcarta L2S transmittances are supplied with first value at
%   the SFC and 101st layer at space. The 101st is above nominal TOA and set to 1.000
%   and is NOT used by Scott's routines. The L2S are swapped for use in the FTC
%   fitting so the 1st value is TOA and 100th value is the SFC (actually below).
%
%   The regression profiles supplied by Klayers in the RTP file are supplied with the
%   first value of prof.ptemp, prof.gas_N etc is at the TOA and the 100th value at SFC.
%
%   Although not used here: SARTA uses a text file for the reference profile that has
%   the first row at the SFC and the 100th row at TOA.
%
%   Before writing out, the data are re-ordered again as follows:
%   roctrans: [nang x nlay x ngas] (6 x 100 x 4) in the case of set 1.
%   gas=1 [1:600], gas=2 [601:1200], gas=3 [1201:1800], gas=4 [1801:2400]
%   thus: for gas=1, ang=1, roctrans(1:6:595, ichan)
%         for gas=1, ang=2, roctrans(2:6:596, ichan)
%         for gas=1, ang=1, roctrans(601:6:1195, ichan) 
%         for gas=1, lev=20, ang=1:6, roctrans(115:120,ichan)
%         etc etc 

% Original created by Scott Hannon, 6 August 1999
% Extensively revised to read in the profile data from RTP; April 2002
% Last updated: 9 May 2002, Scott Hannon - fix bug & check comment
% 13 Apr 2009, S.Hannon - update alltypeband for CrIS; add lfcow2fow
% 30Jan2017 C Hepplewhite: added option to do SAF704 profiles, param cregProf
% 10Mar2018 CLH: updates for prod_2018 (HIT2016): dpath, outd, scang, nang
% 10Dec2018 CLH: prep for prod_2019, new paths.
% 05Feb2020 CLH: Add new sensor CHIRP
% 09Apr2021 CLH: check for 'hi_fcris' or 'fcris' since CHIRP calcs.
% 02Jun2021 CLH: For May2021 kCARTA runs. Update paths. update cris_nsr and nchan
%                assignment is not safe.
% 02Jul2022 CLH: Updated paths for H2020 runs July2022
%                debug old cris_hires checks and myset sets1..5.
%  Oct2022  CLH: option to include all 14 scan angles for fitting the sw (sun)
%                set5. expt to improve CO2 in 4.3 um band.
%  Jan 2024 CLH: Minor changes to accomodate SAF704 profiles. 
%                changed angle subsetting to 8 for thermal sets.
% Jan 2025  CLH: added paths for airs_oco2_pbl modelling
%                added opts.{sensor,prod,build,regset,myset}. for: airs_pbl
%                for nscang=14 (from kCARTA) sets1,2,3 use 7 angles.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/chepplew/projects/sarta/matlabcode
addpath /asl/matlib/h4tools
warning 'off';

%cd /asl/s1/chepplew/projects/sarta/prod_2016/CRIS/Run_trans_cris/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Info for supported typeband
alltypeband={'fowB1', 'fwoB1', 'fmwB2', 'fcowB3', 'fwo_bFsW'};
allngas=    [   4,       4,       3,       5,        4    ]; 
allgasids=zeros( length(alltypeband), max(allngas) );
allgasids(1,1:4)=[2, 3, 1, -2];     % fowB1
allgasids(2,1:4)=[2, 3, 1, -2];     % fwoB1
allgasids(3,1:3)=[2, 6, 1];         % fmwB2
allgasids(4,1:5)=[2, 5, 3, 1, -2];  % fcowB3
allgasids(5,1:4)=[2, 3, 1, -2];     % fwoB3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ok=0;

if (nargin < 3 | nargin > 4)
   error('Unexpected number of input arguements')
end
if (nargin == 3)
   lfcow2fow = 0;
end

% Set nscang (8, 12 or 14 scan angles from L2S OD calcs.)
nscang = opts.nscang;     % 12 or 14;

% Check comment
if (length(comment) > 80)
   error('comment string is too long; max allowed length is 80 char')
   return
end

% hardwire the kCARTA build
build = opts.build; % 'jul2022'; 'may2021'; 'feb2020'; 'mar2018'; 'sep2018'; 'dec2018';

% hardwire the fitting production run
prod_run = ['prod_' opts.prod];   % 2022';

% check mySet
allsets = {'set1','set2','set3','set4','set5'};   % ,'set6','set7'};
myset   = opts.myset;
if(~ismember(myset, allsets)) 
   error('invalid set')
   return
end

% check sensor
csens = upper(opts.csens);
if( ~ismember(csens,{'CRIS_LR','CRIS_HR','AIRS_L1B','AIRS_L1C','IASI', 'CHIRP' ...
     'AIRS_PBL','CRIS_HR_PBL','CHIRP_PBL'}) )
   error('Wrong instrument. Options are: CRIS_{LR,HR},AIRS{ ,L1C,PBL},CHIRP{ ,PBL}')
   return
end

% hardwire number of channels for each sensor
if(strcmp(csens,'CRIS_LR'))  nchan = 1305;  ichan = [1:nchan]; end
if(strcmp(csens,'CRIS_MR'))  nchan = 1683;  ichan = [1:nchan]; end
if(strcmp(csens,'CRIS_HR'))  nchan = 2235;  ichan = [1:nchan]; end
if(strcmp(csens,'AIRS_L1B')) nchan = 2378;  ichan = [1:nchan]; end
if(strcmp(csens,'AIRS_L1C')) nchan = 2834;  ichan = [1:nchan]; end
if(strcmp(csens,'AIRS_PBL')) nchan = 2834;  ichan = [1:nchan]; end
if(strcmp(csens,'IASI'))     nchan = 8461;  ichan = [1:nchan]; end
if(contains(csens,'CHIRP'))    nchan = 1702;  ichan = [1:nchan]; end

% hardwire number of layers in the L2S trans file
nlay = 100;

% check choice of regression profile set
regset = upper(opts.regset);
if( ~ismember(regset,{'R49','SAF704'}) )
  error('Incorrect option for regression profiles. Options are R49, SAF704');
  return
end

% Set the source directories and files:
if(strcmp(regset,'R49'))
  switch build
    case 'jun2016'
      kpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
             'REGR49_400ppm_H2012_June2016/'];
    case 'mar2018'
     %dpath   = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/';  
     dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
            'REGR49_400ppm_H2016_Mar2018/'];
    case 'sep2018'  
      dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
             'REGR49_400ppm_H2016_Sept2018_AIRS2645/'];
    case 'dec2018' 
      dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
             'REGR49_400ppm_H2016_Dec2018_AIRS2834/'];  
    case 'feb2020'
      dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
             'REGR49_400ppm_H2016_Feb2020_AIRS2834_CHIRP/'];
    case 'may2021'
      dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
             'REGR49_400ppm_H2016_May2021_AIRS2834_3CrIS_IASI/'];
    case 'jul2022'
      dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
             'REGR49_400ppm_H2020_July2022_AIRS2834_3CrIS_IASI/'];
      rtpfile   = [dpath 'regr49_1100_400ppm_unitemiss.op.rtp'];
    case 'jan2025a'    % AIRS_OCO2_PBL new LAYER version
      dpath=['/home/sergio/MATLABCODE_Git/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
             'REGR49_400ppm_H2020_Jan2025_PBL_AIRS2834_3CrIS_IASI/'];
      rtpfile   = [dpath 'regr49_pbl.op.rtp'];
  end
  pnums     = [1:48]';
  comment   = [csens ' r49 400ppm H2020 ftc.14a ' build];
  outd_pref = ['/home/chepplew/data/sarta/' prod_run ...
               '/' lower(csens) '/' build '/'];
end
if(strcmp(regset,'SAF704'))
  %dpath   = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/SAF704/';
  dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
         'SAF704_400ppm_H2016_Dec2018_AIRS2834/'];
  rtpfile = [dpath 'save_SAF_704_profiles_29-Apr-2016_1100mb_400ppmv_unitemis.op.rtp'];
  rtpfile = ['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/ECMWF_SAF_137Profiles/' ...
             'save_SAF_704_profiles_29-Apr-2016_1100mb_400ppmv_unitemis.op.rtp'];
  pnums = [1:703]';
  comment = [csens ' SAF704 400ppm CO2 H2016'];
  outd_pref = ['/home/chepplew/data/sarta/' prod_run ...
               '/' lower(csens) '/' build '/'];
end

switch myset
   case 'set1'           % FWOP
     mfile1    = 'F/convolved_kcarta_F_';
     mfile2    = 'FO/convolved_kcarta_FO_';
     mfile3    = 'FWO/convolved_kcarta_FWO_';
     mfile4    = 'FWOP/convolved_kcarta_FWOP_';     % {FWOP, FWOP_CO2_1.0{3,5}; FWOP_Orig}
     mfiles    = {mfile1, mfile2, mfile3, mfile4};
     typeband  = 'fowB1';
     npaths    = 4;
     outd      = [outd_pref 'FWO/'];
     outf      = [csens '_' regset '_allPaths_'];
   case 'set2'           % FOWP
     mfile1    = 'F/convolved_kcarta_F_';
     mfile2    = 'FO/convolved_kcarta_FO_';
     mfile3    = 'FWO/convolved_kcarta_FWO_';
     mfile4    = 'FWOP/convolved_kcarta_FWOP_';
     mfiles    = {mfile1, mfile2, mfile3, mfile4};
     typeband = 'fwoB1';
     npaths    = 4;
     outd      = [outd_pref 'FOW/'];
     outf      = [csens '_' regset '_allPaths_'];
   case 'set3'           % FMW
     mfile1    = 'wvbandF/convolved_kcarta_wvbandF_';   % path 1: (F=1, H2O=0.0, CH4=0.0)
%     mfile2    = 'FO/convolved_kcarta_FO_';             % path 2: (H2O=0.0)
%     mfile3    = 'FWO/convolved_kcarta_FWO_';           % path 3: (all weight = 1)
     mfile2    = 'wvbandFM/convolved_kcarta_FO_';       % path 2: (F=1, H2O=0.0)
     mfile3    = 'wvbandFMW/convolved_kcarta_FWO_';     % path 3: (all weight = 1)
     mfiles    = {mfile1, mfile2, mfile3};
     typeband  = 'fmwB2';
     npaths    = 3;
     outd      = [outd_pref 'wvFMW/'];
%     outd      = [outd_pref 'FMW/'];
     outf      = [csens '_' regset '_allPaths_'];
   case 'set4'
     mfile1    = 'cobandF/convolved_kcarta_cobandF_';   % path 1:
     mfile2    = 'cobandFC/convolved_kcarta_F_';        % path 2:
     mfile3    = 'cobandFCO/convolved_kcarta_FO_';      % path 3:
     mfile4    = 'cobandFCOW/convolved_kcarta_FWO_';    % path 4:
     %mfile5    = 'cobandFCOWP/convolved_kcarta_FWOP_';  % path 5:
     mfile5    = 'FWOP/convolved_kcarta_FWOP_';  % path 5:
     mfiles    = {mfile1, mfile2, mfile3, mfile4, mfile5};
     typeband  = 'fcowB3';
     npaths    = 5;
     outd      = [outd_pref 'FCOW/'];
     outf      = [csens '_' regset '_allPaths_'];
   case 'set5'
     mfile1    = 'F/convolved_kcarta_F_';
     mfile2    = 'FO/convolved_kcarta_FO_';
     mfile3    = 'FWO/convolved_kcarta_FWO_';
     mfile4    = 'FWOP/convolved_kcarta_FWOP_';
     mfiles    = {mfile1, mfile2, mfile3, mfile4};
     typeband  = 'fwo_bFsW';
     npaths    = 4;
     outd      = [outd_pref 'FWOsun/'];
     outf      = [csens '_' regset '_allPaths_'];
   otherwise
     disp('no valid set');
end
      
% hardwire angles (6(8) for sets 1,2,3; 12(14) for sets 4,5,6,7) 
% the kcarta convolved L2S were computed with these 14 angles
secangkc = [[1.00  1.012 1.051 1.19 1.41 1.68 1.99 2.37] ...
            [2.84 3.47 4.30 5.42 6.94 9.02]];
% ORIGINALLY the ftc expects 6 angles for sets 1,2, 3 and 12 for sets 4,5,6,7.
ftcang8  = [1 1.19 1.41 1.68 1.99 2.37 2.84 3.47];
% try 8 annlges from nadir experiment
ftcang8  = [1.00 1.012 1.051 1.19 1.41 1.68 1.99 2.37];
ftcang12 = [[1 1.19 1.41 1.68 1.99 2.37] [2.84 3.47 4.30 5.42 6.94 9.02]];
% To ulitize all fitting angles set1,2,3 use 8 angles and set4,5,6,7 + rem 6.
ftcang14 = [[1 1.012  1.051 1.19 1.41 1.68 1.99 2.37] ...
            [2.84 3.47 4.30 5.42 6.94 9.02]];
% Choose which to use and add id to output convdat file.:
switch nscang
  case 8
    ftcang = ftcang8
    fnosffx = '8a';
  case 12
    ftcang = ftcang12;
    fnosffx = '12a';
  case 14
    ftcang = ftcang14;
    fnosffx = '14a';
end
% subset angles for this combined data - use ysx to subset JA below.
[secang isx ysx] = intersect(ftcang, secangkc);
switch myset
  case {'set1','set2','set3'}
    if(nscang == 12)
      nang = 6;  
      ysx = ysx(1:nang); 
    elseif(nscang == 14)
      nang = 7;
      ysx = ysx(1:nang);
    end
  case {'set4','set5','set6','set7'}
    if(nscang == 8)
      nang = 8;
      ysx = ysx(1:nang);
    elseif(nscang == 12)
      nang = 12;  
      ysx = ysx(1:nang); 
    elseif(nscang == 14)
      nang = 14;
      ysx = ysx(1:nang);
    end
end
secang = ftcang(1:nang);

% print useful summary of what's being processed!
fprintf(1,'Set: %s. Training: %s. Angle set: %d\n', myset,regset,nang);
fprintf(1,'using rtp.ref: %s\n', rtpfile);

% Read the RTP profile data for all 49 profiles
[head, hattr, prof, pattr] = rtpread(rtpfile);       % inherits AIRS header values.
if (head.ptype ~= 1)
   error('rtpfile is not a "layers" profile')
end
% fchan normally comes from the RTP profile but Sergio set the header for AIRS.
%xx    = textread('/home/chepplew/gitLib/ftc_dev/chanLists/list_cris_hrg4');
%fchan = xx(:,2); clear xx;

% Since adding CHIRP the CrIS FSR field names changed:
if(contains(csens, 'CRIS'))
  dtest = load([dpath mfiles{1} int2str(1) '.mat']); 
  if(isfield(dtest,'hi_rcris_all') & contains(csens,'CRIS_HR') )
    ccstr='hi_rcris_all';
    cfstr='hi_fcris'; 
    vchan=dtest.hi_fcris; 
    nchan = length(vchan);
    ichan = [1:nchan];
  end
  if(isfield(dtest,'rcris_all'))
    ccstr='rcris_all';
    cfstr='fcris';
  end
  if(isfield(dtest,'low_rcris_all') & strcmp(csens,'CRIS_LR') ) 
    ccstr='low_rcris_all'; 
    cfstr='low_fcris'; 
    vchan=dtest.low_fcris;
    nchan = length(vchan);
    ichan = [1:nchan];
  end
end

% Loop over the profiles
for ip = [1:length(pnums)]
   clear Y1 Y2 Y3 Y4 J1 J2 J3 J4 JA;
   ipnum=pnums(ip);
   disp(['processing profile ' int2str(ipnum)])

   %%%%%%%%%%%%%%%%%%%%%%%%%
   % Read in the matlab data
   % The matfile must contain:
   %       ctrans(nchan x [nang,nlay,nsets]) {layers bot to top}
   %       secang(nang)
   %       fchan(nchan)
   %       ichan(nchan)
   %       typeband
   %       nchan, nang, nlay, nsets

%   djunk=dir(strcat(dpath1,mfile));
%   [irow,icol]=size(djunk);
%   if (irow ~= 1)
%      error(['Error openning matlab file ' mfile])
%   end
%   clear djunk irow icol
% collect the data and re-arrange before concatenating
   if(numel(mfiles) >=1 )
     Y1 = load([dpath mfiles{1} int2str(ipnum) '.mat']); 
     if(contains(csens,'CRIS'));  J1 = permute(Y1.(ccstr),[2 1 3]); end; 
     if(contains(csens,'AIRS'));  J1 = permute(Y1.rairs_all,[2 1 3]); end;
     if(contains(csens,'IASI'));  J1 = permute(Y1.riasi_all,[2 1 3]); end;
     if(contains(csens,'CHIRP')); J1 = permute(Y1.med_rcris_all,[2 1 3]); end;
   end;
   if(numel(mfiles) >=2 )
     Y2 = load([dpath mfiles{2} int2str(ipnum) '.mat']); 
     if(contains(csens,'CRIS'));  J2 = permute(Y2.(ccstr),[2 1 3]); end;
     if(contains(csens,'AIRS'));  J2 = permute(Y2.rairs_all,[2 1 3]); end;
     if(contains(csens,'IASI'));  J2 = permute(Y2.riasi_all,[2 1 3]); end;
     if(contains(csens,'CHIRP')); J2 = permute(Y2.med_rcris_all,[2 1 3]); end;
   end;
   if(numel(mfiles) >=3 )
     Y3 = load([dpath mfiles{3} int2str(ipnum) '.mat']); 
     if(contains(csens,'CRIS'));  J3 = permute(Y3.(ccstr),[2 1 3]); end;
     if(contains(csens,'AIRS'));  J3 = permute(Y3.rairs_all,[2 1 3]); end;
     if(contains(csens,'IASI'));  J3 = permute(Y3.riasi_all,[2 1 3]); end;
     if(contains(csens,'CHIRP')); J3 = permute(Y3.med_rcris_all,[2 1 3]); end;
   end;
   if(numel(mfiles) >= 4)
     Y4 = load([dpath mfiles{4} int2str(ipnum) '.mat']); 
     if(contains(csens,'CRIS'));  J4 = permute(Y4.(ccstr),[2 1 3]); end;
     if(contains(csens,'AIRS'));  J4 = permute(Y4.rairs_all,[2 1 3]); end;
     if(contains(csens,'IASI'));  J4 = permute(Y4.riasi_all,[2 1 3]); end;
     if(contains(csens,'CHIRP')); J4 = permute(Y4.med_rcris_all,[2 1 3]); end;
   end;
   if(numel(mfiles) >=5)
     Y5 = load([dpath mfiles{5} int2str(ipnum) '.mat']); 
     if(contains(csens,'CRIS'));  J5 = permute(Y5.(ccstr),[2 1 3]); end;
     if(contains(csens,'AIRS'));  J5 = permute(Y5.rairs_all,[2 1 3]); end;
     if(contains(csens,'IASI'));  J5 = permute(Y5.riasi_all,[2 1 3]); end;
     if(contains(csens,'CHIRP')); J5 = permute(Y5.med_rcris_all,[2 1 3]); end;
   end;
   if(numel(mfiles) >=6)
     fprintf(1,'ERROR: cant deal with more than 4 mixed paths\n'); end

   % Re-arrange data and concatenate.
   if(npaths == 3) JA = cat(4, J1, J2, J3); end
   if(npaths == 4) JA = cat(4, J1, J2, J3, J4); end
   if(npaths == 5) JA = cat(4, J1, J2, J3, J4, J5); end

%  size(JA): nchans x  12(14) x 101 x 4
%  take first 6(8) angles as required by fitftc for sets 1 to 3
   if(ismember(myset,{'set1','set2','set3'}))
     disp([myset ': retain ' num2str(nang) '  angles']);
     JA = JA(:, ysx, :, :);
   end
   if(ismember(myset,{'set4','set5','set6','set7'}))
     disp([myset ': retain ' num2str(nang) '  angles']);
     JA = JA(:, ysx, :, :);
   end
% for Scott's fitting code must have 100 layers, Sergio's 101th layer is set to 1.
   if(size(JA,3) == 101) 
     fprintf(1,'Original L2S files have 101 layers: truncate to 100 for fitftc\n');
     JA = JA(:,:,1:100,:);
   end
   if(contains(csens,'CRIS'));  fchan = Y1.(cfstr); end
   if(contains(csens,'AIRS'));  fchan = Y1.fairs; end
   if(contains(csens,'IASI'));  fchan = Y1.fiasi; end
   if(contains(csens,'CHIRP')); fchan = Y1.med_fcris; end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Check typeband and determine gas IDs
   ii = strmatch(typeband,alltypeband);
   if ( length(ii) > 1 | ii < 1)
      ii
      typeband
      alltypeband
      ipnum
      error(['invalid type_band: ' typeband]);
   end
   if (ii ~= 4)
      lfcow2fow = 0;
   end
   ngas=allngas(ii);
   if (ngas ~= npaths)
     ngas
     nsets
     ipnum
     error('mismatch in expected ngas and matfile nsets')
   end
   gasids=allgasids(ii,1:ngas);

%{
figure(2);clf;h1=semilogy(squeeze(Y3.rcris_all(1,1191,:)),prof.plevs(1:end-1,1),'.-');grid on;
  ax=gca; ax.YDir='reverse';ylim([1 1030]);
figure(2);clf;hold on; for i=1:19:100 plot(fchan,JA(:,1,i,1),'-');end; grid on;
figure(2);clf;hold on; for i=1:4 plot(squeeze(JA(714,1,:,i)),[1:100],'.-');end; grid on;
figure(2);clf;semilogy(squeeze(JA(1312,1,:,3)),prof.plevs(1:end-1,1),'.-');grid on;
  ax=gca; ax.YDir='reverse';ylim([1 1030]);xlabel('FWO rcris.all');ylabel('pressure hPa');
  title('convolved.kcarta.FMW.1.mat, sec=1, 1576wn');
%}
   %%%%%%%%%%%%%%%%%%
   % the Klayers data are in correct layer order (1st=TOA) 
   % Load temperature
   ii=prof.nlevs(ipnum) - 1; % number of layers
   if (ii ~= nlay)
      ii
      nlay
      ipnum
      error('mismatch in number of layers')
   end
   temp=prof.ptemp(1:nlay,ipnum);


   %%%%%%%%%%%%%%%%%%
   % Load gas amounts
   amount=zeros(nlay,ngas);
   for ig=1:ngas

      ii=intersect( head.glist, abs(gasids(ig)) );
      if (length(ii) ~= 1)
         ipnum
         gasids(ig)
         head.glist
         error('required gas not found in RTP file')
      end
      ii=abs(gasids(ig));
      eval(['amount(:,ig)=prof.gas_' int2str(ii) '(1:nlay,ipnum);'])
   end
   % Convert amount from molecules/cm^2 to kmoles/cm^2
   amount=amount/6.02214199E+26;
%{
figure(4);clf; loglog(amount(:,1),prof.plevs(1:100,1),'.-');  
figure(4);clf; loglog(temp,prof.plevs(1:100,1),'.-');
%}
   % Assign dummy res and rnfwhm
   res      = zeros(nchan,1);
   rnfwhm   = zeros(nchan,1);
   roctrans = [];

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Re-order the ctrans layers to run top to bottom
   if (ip ==1)
      indr       = nlay:-1:1;
      nanglaygas = nang*nlay*ngas;
      nanglay    = nang*nlay;
      indis      = 1:nanglaygas;
      indwant    = indis;
      ilay   = 1:nlay;
      ileft  = nang*(ilay-1);
      iright = nang*(nlay - ilay);
      for igas=1:ngas
         ioffset = (igas - 1)*nanglay;
         for ia=1:nang
            indwant(ioffset + ileft + ia)=indis(ioffset + iright + ia);
         end
      end
   end
   roctrans=JA(:,indwant)';       % [2400x2235] w/6 angs 4 paths. [6000x2235] w/12 angles 5 paths
%{
figure(4);clf;hold on; for i=115:120 plot(fchan,roctrans(i,:),'-');end;grid on; %near TOA
figure(5);clf;hold on; for i=415:420 plot(fchan,roctrans(i,:),'-');end;grid on; %
figure(6);clf;hold on; plot(roctrans(1:6:595,54),[1:100],'-');grid on;
%}

fprintf(1,'lfcow2fow= %5i  ngas= %5i\n', lfcow2fow,ngas);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Write out the merged profile + transmittance data file
   if(exist(outd,'dir') ~= 7) 
     disp([outd ' does not exist, creating it']); 
     ires = mkdir(outd);
     if(ires == 0) error('could not create directory'); end
   end

   outname=[outd outf fnosffx '_' int2str(ipnum) '.dat'];
   titlecom=[comment ' p#' int2str(ipnum)];
   disp(['Writing data to file: ' outname]);
   
   if (lfcow2fow == 1 )
     [iok]=wrt_convdat_fcow2fow(outname, nang, nlay, ngas, nchan, gasids, secang,...
     titlecom, temp, amount, fchan, ichan, res, rnfwhm, roctrans);
   else
      [iok]=wrt_convdat(outname, nang, nlay, ngas, nchan, gasids, secang, ...
        titlecom, temp, amount, fchan, ichan, res, rnfwhm, roctrans);
   end
   if (iok == 0)
      error(['Error detected writing ' outname]);
   end

end
ok=1;

%%% end of function %%%
