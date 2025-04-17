% make_r49_2sfc_8view_2sol_test.m
% 
% Take original 49 training profiles and make test set with:
%   2 surfaces (black + sea), 8 view angles, solar noon and night.
%   49 x 2 x 8 x 2 = 1568
% 

cd /home/chepplew/projects/sarta/matlabcode

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools/
addpath /home/chepplew/myLib/matlib/rtptools1   % replicate_rtp_headprof
addpath /asl/matlib/aslutil                     % int2bits
addpath /home/chepplew/gitLib/rtp_prod2/emis    % emis_sea

% Which LAYERING to use (airs.plevs or pbl.plevs)
opts.plevs = 'pbl';

% klayers execs
klayersexe.pbl = '/home/chepplew/gitLib/klayersV205/BinV201/klayers_pbl_wetwater_test';
klayersexe.air = '/asl/packages/klayersV205/BinV201/klayers_airs_wetwater';

% assign view and solar angles to use
satzen = [0.0 8.8322 17.9223 32.8244 44.8285  53.4704  59.8336  65.0428];
solzen = [0.0 150.0]; 

% load sea emissivity at nadir, no wind.
wspeed = 0.0;
zang   = 0.0;
[nemis, efreq, seaemis] = emis_sea(zang, wspeed);

fn.ip_rtp = '/home/chepplew/projects/klayers_wrk/regr49_1100.ip.rtp';

% load basic r49 profiles LAYERS file
switch opts.plevs
  case 'pbl'
    fn.r49_rtp = '/home/chepplew/data/sarta/prod_2025/generic/regr49_pbl.op.rtp';
  case 'airs'
    fn.r49_rtp = ['/home/sergio/MATLABCODE_Git/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
    'REGR49_400ppm_H2020_July2022_AIRS2834_3CrIS_IASI/regr49_1100_400ppm_unitemiss.op.rtp'];
end

[head,hatt,prof,patt] = rtpread(fn.ip_rtp);
% set CO2 surface abundance to 610 ppm mass mixing ratio (400 ppmv).
co2_rat = 610.0E-6/max(prof.gas_2(1:9,1));
prof.gas_2 = co2_rat*prof.gas_2;

% Three stage replication for view angles, surfaces and solar.
% step 1: replicate for 8 view angles
h2 = struct;
p2 = struct;
for ii = 1 : 49
   [hy,py] = replicate_rtp_headprof(head,prof,ii,8);
   py.emis(:,1:8) = repmat(ones(19,1),1,8);
   py.rho(:,1:8)  = repmat(zeros(19,1),1,8);
   for jj = 1:8
     py.satzen(jj) = satzen(jj);
   end
   if ii == 1
     h2 = hy;
     p2 = py;
   else
     [h2,p2] = cat_rtp(h2,p2,hy,py);
   end
   fprintf(1,'.')
end
fprintf(1,'\n')
clear hy py;

% step 2: replicate again for 2 surfaces
% 1:2:end = black. 2:2:end = sea.
h3 = struct;
p3 = struct;
for ii = 1:length(p2.satzen)
   [hy,py] = replicate_rtp_headprof(h2,p2,ii,2);
   py.emis(:,2) =  seaemis;
   py.rho(:,2)  = (1 - seaemis)/pi;
   if ii == 1
     h3 = hy;
     p3 = py;
   else
     [h3,p3] = cat_rtp(h3,p3,hy,py);
   end
   if(~mod(ii,10)) fprintf(1,'.'); end
end
fprintf(1,'\n')
clear hy py;

% step 3: replicate again for 2 solar (noon, night)
% 1:2:end = night. 2:2:end = noon. 
h4 = struct;
p4 = struct;
for ii = 1:length(p3.satzen)
   [hy,py] = replicate_rtp_headprof(h3,p3,ii,2);
   py.solzen(:,2) = 0.0;
   if ii == 1
     h4 = hy;
     p4 = py;
   else
     [h4,p4] = cat_rtp(h4,p4,hy,py);
   end
   if(~mod(ii,20)) fprintf(1,'.'); end
end
fprintf(1,'\n')
clear hy py;

% summary:
% all unique satzens are at indexes: 1*8*2*2:end = X:32:end, X=[1,5,9,13,17,21,25,29];
% all nadir views at noon over black surface:
%   iix = find(p4.satzen==0 & p4.solzen==0 & p4.emis(1,:) == 1);
% all 8 view angles at noon over the sea:
%   iix = find(p4.solzen == 0 & p4.emis(1,:) < 0.9999);
% r49 at nadir & night & black.
%   iix = [1:32:nprofs];

% boundary conditions
h4.pmax   = 1045.0;
p4.spres  = 1013.25*ones(size(p4.spres));
p4.salti  = zeros(size(p4.salti));
p4.co2ppm = 400.0*ones(size(p4.co2ppm));

% ==== prep layers file for AIRS.L1c 2834 channel set =======
ichan = hdfread('/home/chepplew/myLib/data/airs_l1c_srf_tables_lls_new.hdf','chanid');
vchan = hdfread('/home/chepplew/myLib/data/airs_l1c_srf_tables_lls_new.hdf','freq');
% load('/home/chepplew/myLib/data/fcris_hires_4grd.mat','ichan');;
% load('/home/chepplew/myLib/data/fcris_hires_4grd.mat','vchan');
% -------------------------------------------------------
% Get LAYERS using klayers w/ PBL layering
tempPath = mktemp();
fn.rtp_in = [tempPath '_test_profs_ip.rtp'];
rtpwrite(fn.rtp_in,h4,hatt,p4,patt);
%
fn.rtp_op1 = [tempPath '_test_profs_pbl_op.rtp'];

command1 = [klayersexe.pbl ' fin=' fn.rtp_in ' fout=' fn.rtp_op1  ...
      ' > /home/chepplew/data/scratch/kla_out.log'];
tic
system(command1)
toc

[hda, haa, pda,paa] = rtpread(fn.rtp_op);
hda.vchan = single(vchan(:));
hda.ichan = int32(ichan(:));
hda.nchan = length(ichan);
% Save RTP file.
sav.dir = '/home/chepplew/data/sarta/prod_2025/generic/';
sav.fn  = 'r49_1013_400ppm_8angs_2sfc_2sols_pbl_2834_op.rtp'
rtpwrite([sav.dir sav.fn], hda, hatt, pda, patt);
% -------------------------------------------------------
% Get LAYERS using klayers w/ PBL layering
fn.rtp_op2 = [tempPath '_test_profs_airslay_op.rtp'];
command2 = [klayersexe.air ' fin=' fn.rtp_in ' fout=' fn.rtp_op2  ...
      ' > /home/chepplew/data/scratch/kla_out.log'];
tic
system(command2)
toc

[hda, haa, pda,paa] = rtpread(fn.rtp_op);
hda.vchan = single(vchan(:));
hda.ichan = int32(ichan(:));
hda.nchan = length(ichan);
% Save RTP file.
sav.dir = '/home/chepplew/data/sarta/prod_2025/generic/';
sav.fn  = 'r49_1013_400ppm_8angs_2sfc_2sols_airslay_2834_op.rtp'
rtpwrite([sav.dir sav.fn], hda, hatt, pda, patt);
% ----------------------------------------------------------

% Check layers file
[hda, haa, pda,paa] = rtpread(fn.rtp_op1);   % fn.rtp_op2
nprofs = length(pda.spres);

clear nlevs
for i=1:nprofs
  nlevs(i) = pda.nlevs(i)-1;
end
% r49 at nadir & night & black.
iix=[1:32:nprofs];
figure(2);clf
i=1;
  semilogy(pda.ptemp([1:nlevs(i)],iix(i)), pda.plevs([1:nlevs(i)],iix(i)),'.-');
  hold on;
for i=2:length(iix)
  semilogy(pda.ptemp([1:nlevs(i)],iix(i)), pda.plevs([1:nlevs(i)],iix(i)),'.-');
  hold on;
end

% IF satisfied - save profiles to rtp file
disp(['saving test profiles to: ' [sav.dir sav.fn]]);
rtpwrite([sav.dir sav.fn], h4,hatt,p4,patt);

% END




%{
% This matrix gives the order of profiles after duplication and perturbation
% at nadir ONLY:
  pind =  {{'WV',   1,  [1:9:433]},...
           {'CO2',  2,  [2:9:434]},...
           {'O3',   3,  [3:9:435]},...
           {'N2O',  4,  [4:9:436]},...
           {'CO',   5,  [5:9:437]},...
           {'CO',   6,  [6:9:438]},...
           {'SO2',  9,  [7:9:439]},...
           {'HN3', 11,  [8:9:440]},...
           {'NHO3',12,  [9:9:441]},...
           {'unp',  0,  [442:490]}};

% Define output directory:
outdr = '/home/chepplew/data/sarta/prod_2019/generic/';

% original regr49 profiles:
srcdr  = ['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
          'REGR49_400ppm_H2016_Dec2018_AIRS2834/'];
srcrtp = [srcdr 'regr49_1013_400ppm_unitemiss.op.rtp'];

fnrtp =['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
        'REGR49_400ppm_H2016_Mar2018_NH3/stdNH3_1100mb_op_400ppm.rtp'];

[head, ~, prof, ~]       = rtpread(srcrtp);
[htmp, ~, ptmp, ~]       = rtpread(fnrtp);

% some fields are missing from head

head.glist = htmp.glist;
head.gunit = htmp.gunit;
head.ngas  = htmp.ngas;

% Do the nadir duplication
h2 = struct;
p2 = struct;
for ii = 1 : 49
   [hy,py] = replicate_rtp_headprof(head,prof,ii,4);
   py.emis(:,1) = ones(19,1);

   if ii == 1
     h2 = hy;
     p2 = py;
   else
     [h2,p2] = cat_rtp(h2,p2,hy,py);
   end
   fprintf(1,'.')
end
fprintf(1,'\n')

% append the unperturbed original profile set to the end
% complicated 'cos only subset of fields were duplicated
fnames2 = fieldnames(p2);
fnames1 = fieldnames(prof);
diff(fnames1, fnames2);

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

% Force all profiles to have unit surface emissivity
prof2.emis = ones(size(prof2.emis));

% write out pertubation set:
outdr = '/home/chepplew/data/sarta/prod_2019/generic/';
fortp = 'r49_1100_400p_seaemis_nadir_pert_v2.rtp';

rtpwrite([outdr fortp], head2, hatt, prof2, patt);

% -------------------------------------
%          Duplicate for non-nadir
% -------------------------------------
angles = [0  8.8322 17.9223 32.8244 44.8285  53.4704  59.8336  65.0428];

h3 = struct;
p3 = struct;
nprof = size(prof2.satzen,2);     % expecting 490
for ii = 1 : nprof
   [hy,py] = replicate_rtp_headprof(head2,prof2,ii,8);
   for jj = 2:8
     py.satzen(jj) = angles(jj);
   end
  if ii == 1
    h3 = hy;
    p3 = py;
  else
    [h3, p3] = cat_rtp(h3, p3, hy, py);
  end
  if(~mod(ii,49)); fprintf(1,'.'); end
end

% retain same prof.emis = 1, and head structure

% save rtp file
fortp = 'r49_1100_400p_seaemis_8angs_pert_v2.rtp';

rtpwrite([outdr fortp], h3, hatt, p3, patt);

% ----------------------------------------
% Algorithm for indexing each perturbation
% ----------------------------------------

indx = struct;
for ii=1:49 
  for kk=1:8 
    jj=8*(ii-1)+kk;
    indx.wv(jj)   = kk    + (ii-1)*72;
    indx.co2(jj)  = 8+kk  + (ii-1)*72;
    indx.o3(jj)   = 16+kk + (ii-1)*72;
    indx.n2o(jj)  = 24+kk + (ii-1)*72;
    indx.co(jj)   = 32+kk + (ii-1)*72;
%    indx.cox(jj)  = 32+kk + (ii-1)*72;
    indx.ch4(jj)  = 40+kk + (ii-1)*72;
    indx.so2(jj)  = 48+kk + (ii-1)*72;
    indx.nh3(jj)  = 56+kk + (ii-1)*72;
    indx.hno3(jj) = 64+kk + (ii-1)*72;
    %disp([num2str(jj) ': ' num2str(indx.so2(jj))]);
    
  end 
end
indx.unp = [3529:3920];



%--------------------- END -------------------

% Check indexing for select gas at select angle
% 1.1: unperturbed at nadir:
izo = indx.unp(1:8:392);
% 1.2 perturbed water at nadir:
xzo = indx.wv(1:8:392);
figure(11);clf;plot(p3.satzen(xzo),'o');
iip = [1,2,49];
figure(11);clf;
 plot(p3.gas_1(:,[xzo(iip),izo(iip)]), p3.plevs(:,[xzo(iip),izo(iip)]),'-');
 grid on; set(gca,'YDir','reverse');title('H2O');


% H2O
ipp = pind{1}{3};
ipu = pind{10}{3};
ip = [1,2,49];
figure(1);clf;
 plot(prof2.gas_1(:,[ipp(ip),ipu(ip)]), prof2.plevs(:,[ipp(ip),ipu(ip)]),'-');
 grid on; set(gca,'YDir','reverse');title('H2O');

% CO2
ipp = pind{2}{3};
ipu = pind{10}{3};
ip = 1;
figure(1);clf;
 plot(prof2.gas_2(:,[ipp(ip),ipu(ip)]), prof2.plevs(:,[ipp(ip),ipu(ip)]),'-');
 grid on; set(gca,'YDir','reverse');title('CO2');

% O3
ipp = pind{3}{3};
ipu = pind{10}{3};
ip = 1;
figure(1);clf;
 plot(prof2.gas_3(:,[ipp(ip),ipu(ip)]), prof2.plevs(:,[ipp(ip),ipu(ip)]),'-');
 grid on; set(gca,'YDir','reverse');title('O3');

% N2O
ipp = pind{4}{3};
ipu = pind{10}{3};
ip = 1;
figure(1);clf;
 plot(prof2.gas_4(:,[ipp(ip),ipu(ip)]), prof2.plevs(:,[ipp(ip),ipu(ip)]),'-');
 grid on; set(gca,'YDir','reverse');title('N2O');

% CO
ipp = pind{5}{3};
ipu = pind{10}{3};
ip = 1;
figure(1);clf;
 plot(prof2.gas_5(:,[ipp(ip),ipu(ip)]), prof2.plevs(:,[ipp(ip),ipu(ip)]),'-');
 grid on; set(gca,'YDir','reverse');title('CO')

% CH4
ipp = pind{6}{3};
ipu = pind{10}{3};
ip = 1;
figure(1);clf;
 plot(prof2.gas_6(:,[ipp(ip),ipu(ip)]), prof2.plevs(:,[ipp(ip),ipu(ip)]),'-');
 grid on; set(gca,'YDir','reverse');title('CH4');

% SO2
ipp = pind{7}{3};
ipu = pind{10}{3};
ip = 1;
figure(1);clf;
 plot(prof2.gas_9(:,[ipp(ip),ipu(ip)]), prof2.plevs(:,[ipp(ip),ipu(ip)]),'-');
 grid on; set(gca,'YDir','reverse');title('SO2');

% NH3
ipp = pind{8}{3};
ipu = pind{10}{3};
ip = 1;
figure(1);clf;
 plot(prof2.gas_11(:,[ipp(ip),ipu(ip)]), prof2.plevs(:,[ipp(ip),ipu(ip)]),'-');
 grid on; set(gca,'YDir','reverse');title('NH3');

% HNO3
ipp = pind{9}{3};
ipu = pind{10}{3};
ip = 1;
figure(1);clf;
 plot(prof2.gas_12(:,[ipp(ip),ipu(ip)]), prof2.plevs(:,[ipp(ip),ipu(ip)]),'-');
 grid on; set(gca,'YDir','reverse');title('HNO3');


%}
