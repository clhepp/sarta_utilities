% perturb49_kc_vs_sarta.m

% comparison of TOA radiances derived from kcarta and SARTA for hi-res CrIS
% with minor gas profiles perturbed as follows:
% profiles:
%
% rev.1 Simple 10% column perturbations all gases
% rev 2: from LLS slack:sarta_upgrade 1Mar2019 "CO2: 10 ppm,  CH4: 3% , 
% N2O: 3%  CO, 10% is good, run at 30% too so can see how errors vary.  
% HNO3: not sure yet,  HDO: try 50% as stated"
%
% 

cd /home/chepplew/projects/sarta/matlabcode

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools/
addpath /asl/matlib/aslutil            % int2bits
addpath /home/chepplew/myLib/matlib/rtptools      % replcate_rtp_headprof

% hardwire version control
vers = 1;

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
outdr = '/home/chepplew/data/sarta/prod_2020/generic/';

% original regr49 profiles:
srcdr  = ['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
          'REGR49_400ppm_H2016_Dec2018_AIRS2834/'];
srcrtp = [srcdr 'regr49_1013_400ppm_unitemiss.op.rtp'];

srcrtp =['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
        'REGR49_400ppm_H2016_Mar2018_NH3/stdNH3_1100mb_op_400ppm.rtp'];

srcrtp = ['/home/chepplew/data/sarta/prod_2025/generic/regr49_1100_pbl_nh3.op.rtp'];


[head, hatt, prof, patt] = rtpread(srcrtp);
%[htmp, hatt, prof, patt] = rtpread(fnrtp);
%{
% some fields are missing from head
head.glist = htmp.glist;
head.gunit = htmp.gunit;
head.ngas  = htmp.ngas;
%}
switch vers
  case 1
% *************** Do the nadir perturbation vers.1 *****************
h2 = struct;
p2 = struct;
for ii = 1 : 49
   [hy,py] = replicate_rtp_headprof(head,prof,ii,9);
   py.gas_1(:,1)  = 1.1   * prof.gas_1(:,ii);          % H2O
   py.gas_2(:,2)  = 1.05 * prof.gas_2(:,ii);           % CO2 20ppm = 5%
   py.gas_3(:,3)  = 1.1   * prof.gas_3(:,ii);          % O3
   py.gas_4(:,4)  = 1.05  * prof.gas_4(:,ii);          % N2O 5%
   py.gas_5(:,5)  = 1.1   * prof.gas_5(:,ii);          % CO 10%
   py.gas_6(:,6)  = 1.1   * prof.gas_6(:,ii);          % CH4
   py.gas_9(:,7)  = 1.1   * prof.gas_9(:,ii);          % SO2
   py.gas_11(:,8) = 1.1   * prof.gas_11(:,ii);         % NH3
   py.gas_12(:,9) = 1.1   * prof.gas_12(:,ii);         % HNO3
   if ii == 1
     h2 = hy;
     p2 = py;
   else
     [h2,p2] = cat_rtp(h2,p2,hy,py);
   end
   fprintf(1,'.')
end
fprintf(1,'\n')

  case 2
% *************  Do the nadir perturbation vers.2 *******************
h2 = struct;
p2 = struct;
for ii = 1 : 49
   [hy,py] = replicate_rtp_headprof(head,prof,ii,9);
   py.gas_1(:,1)  = 1.1   * prof.gas_1(:,ii);
   py.gas_2(:,2)  = 1.025 * prof.gas_2(:,ii);          % CO2 10ppm = 2.5%
   py.gas_3(:,3)  = 1.1   * prof.gas_3(:,ii);
   py.gas_4(:,4)  = 1.03  * prof.gas_4(:,ii);          % N2O 3%
   py.gas_5(:,5)  = 1.1   * prof.gas_5(:,ii);          % CO 10%
   py.gas_5(:,6)  = 1.3   * prof.gas_5(:,ii);          % CO 30%
   py.gas_6(:,7)  = 1.03  * prof.gas_6(:,ii);          % CH4 3%
   py.gas_9(:,8)  = 1.1   * prof.gas_9(:,ii);
   py.gas_11(:,9) = 1.1   * prof.gas_11(:,ii);
%   py.gas_12(:,9) = 1.1  * prof.gas_12(:,ii);
   if ii == 1
     h2 = hy;
     p2 = py;
   else
     [h2,p2] = cat_rtp(h2,p2,hy,py);
   end
   fprintf(1,'.')
end
fprintf(1,'\n')

end    %  END switch version

% append the unperturbed original profile set to the end
% complicated 'cos only subset of fields were duplicated
fnames2 = fieldnames(p2);
fnames1 = fieldnames(prof);
setdiff(fnames1, fnames2);

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

% write out pertubation set for NADIR:
outdr = '/home/chepplew/data/sarta/prod_2019/generic/';
outdr = '/home/chepplew/data/sarta/prod_2025/generic/';
fortp = ['r49_1100_400p_seaemis_nadir_pert_v' num2str(vers) '.rtp'];

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
fortp = ['r49_1100_400p_unitemis_8angs_gas_pert_v' num2str(vers) '.rtp'];
fortp = ['r49_1013_400p_pbl_unitemis_8angs_gas_pert_v' num2str(vers) '.rtp'];

rtpwrite([outdr fortp], h3, hatt, p3, patt);

% ----------------------------------------
% Algorithm for indexing each perturbation
% ----------------------------------------

switch vers
  case 1
indx = struct;
for ii=1:49 
  for kk=1:8 
    jj=8*(ii-1)+kk;
    indx.wv(jj)   = kk    + (ii-1)*72;
    indx.co2(jj)  = 8+kk  + (ii-1)*72;
    indx.o3(jj)   = 16+kk + (ii-1)*72;
    indx.n2o(jj)  = 24+kk + (ii-1)*72;
    indx.co(jj)   = 32+kk + (ii-1)*72;
    indx.ch4(jj)  = 40+kk + (ii-1)*72;
    indx.so2(jj)  = 48+kk + (ii-1)*72;
    indx.nh3(jj)  = 56+kk + (ii-1)*72;
    indx.hno3(jj) = 64+kk + (ii-1)*72;
    %disp([num2str(jj) ': ' num2str(indx.so2(jj))]);
    
  end 
end
indx.unp = [3529:3920];

  case 2
  
end   % END switch vers

%--------------------- END -------------------

%{
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
