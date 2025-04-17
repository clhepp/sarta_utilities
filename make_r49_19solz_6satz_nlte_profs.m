% make_r49_19solz_6satz_nlte_profs.m

% Make test set for nonLTE using regr49 set
% Based on make_r49_perturb_profs.m
%
% Option: airslay or PBL layers
% Jan 2025.  CLH; Partial sets of options currently implemented.

cd /home/chepplew/projects/sarta/matlabcode

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools/
addpath /asl/matlib/aslutil            % int2bits
addpath /home/chepplew/myLib/matlib/rtptools      % replcate_rtp_headprof

% check options
if(~all(isfield(opts,{'csens','prod','build','lays','regset','vers'})))
  error('please include all required options')
  return;
end

% version string
vers = lower(opts.vers);

% AIRSLAY or PBL Layers
if(~ismember(opts.lays,{'airs','pbl'}))
  error('invalid layer option')
  return;
end
lays = lower(opts.lays);

% solar zenith angles to use:
solz = [0:10:80 85 87 90 92 94 96 98 100 105 120];
nsolz = length(solz);

% Satellite zenith angles to use:
satz = [0:10:60];
nsatz = length(satz);

% flat sea surface emissivity to use
% TBD

% original regr49 profiles:
srcdr  = ['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
          'REGR49_400ppm_H2016_Dec2018_AIRS2834/'];

switch lower(opts.lays)
  case 'pbl'
    srcrtp = ['/home/chepplew/data/sarta/prod_2025/generic/regr49_1100_pbl_nh3.op.rtp'];
  case 'airs'
    srcrtp =['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
             'REGR49_400ppm_H2016_Mar2018_NH3/stdNH3_1100mb_op_400ppm.rtp'];
  otherwise
    error('cannot process lays option')
    return
end

% load base set of profiles (regr49)
[head, hatt, prof, patt] = rtpread(srcrtp);
nprofs = length(prof.solzen);
%[htmp, hatt, prof, patt] = rtpread(fnrtp);
%{
% some fields are missing from head
head.glist = htmp.glist;
head.gunit = htmp.gunit;
head.ngas  = htmp.ngas;
%}

% *************  Do the nadir solar zenith var *******************
h1 = struct;
p1 = struct;
for ii = 1 : nprofs     % 49
   [hy,py] = replicate_rtp_headprof(head,prof,ii,19);
   for jj = 1:nsolz
      py.solzen(jj)  = solz(jj);
   end
   if ii == 1
     h1 = hy;
     p1 = py;
   else
     [h1,p1] = cat_rtp(h1,p1,hy,py);
   end
   fprintf(1,'.')
end
fprintf(1,'\n')
clear hy py;

% ------------ Do the satellite zenith var --------------
nprofs = length(p1.solzen);
h2 = struct;
p2 = struct;
for ii = 1: nprofs
  [hy, py] = replicate_rtp_headprof(h1,p1,ii,nsatz);
  for jj = 1:nsatz
     py.satzen(jj) = satz(jj);
  end
  if ii == 1
    h2 = hy;
    p2 = py;
  else
    [h2, p2] = cat_rtp(h2, p2, hy, py);
  end
  if(~mod(ii,20)) fprintf(1,'.'); end
end
fprintf(1,'\n')
clear hy py;

%{
% template replicator
ho = struct;
po = struct;
for ii = 1:49
  [hx, px] = replicate_rtp_headprof(head,prof,ii,1);
  if ii == 1
    ho = hx;
    po = px;
  else
    [ho, po] = cat_rtp(ho, po, hx, px);
  end
end

head2 = struct;
prof2 = struct;
[head2,prof2] = cat_rtp(h2,p2, ho, po);
%}

% Force all profiles to have unit surface emissivity
p2.emis = ones(size(p2.emis));

% write out pertubation set for NADIR:
outdr = ['/home/chepplew/data/sarta/' prod_run '/generic/'];
switch opts.lays
  case 'airs'
    fortp = ['r49_1100_400p_unitemis_19solz_6satz_' num2str(opts.vers) '.rtp'];
  case 'pbl'
    fortp = ['r49_1013_400p_pbl_unitemis_19solz_6satz_' num2str(opts.vers) '.rtp'];
  otherwise
    error('can not process opts.lays')
    return
end
disp(['Saving rtp file to: ' [outdr fortp]])

rtpwrite([outdr fortp], h2, hatt, p2, patt);

% ----------------------------------------
% Subsets for indexing each solzen v satzen
% ----------------------------------------

clear indx;
  for jj = 1:19
    for kk = 1:6 
      indx{jj}{kk} = find(p2.solzen == solz(jj) & p2.satzen == satz(kk));
    end
  end


%--------------------- END -------------------

