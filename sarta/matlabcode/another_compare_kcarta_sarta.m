cd /home/chepplew/projects/sarta/prod_2018/

addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil                       % rad2bt

myvers   = 'cris_hr';
csarta   = 'crisg4_may18_nh3';
ufo      = ['./sarta_data/' myvers '/sar_' csarta '_r49_unpert_seaemis.rtp'];
SARTAEXE = '/home/chepplew/gitLib/sarta/bin/crisg4_may18_nh3';

rtp49 = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP',...
         '/stdNH3_1100mb_op_400ppm.rtp'];
rtp49 = ['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2016_Mar2018/',...
          'regr49_1013_400ppm_unitemiss.op.rtp'];
[hd ha pd pa] = rtpread(rtp49);

load('/home/chepplew/projects/sarta/prod_2016/cris_hr/rtp_head_str_2235g4.mat');
hdr2.pfields = 1;
hdr2.glist = sort([hdr2.glist; 11]);
hdr2.gunit = [hdr2.gunit; 1];

rtp49x = ['./sarta_data/' myvers '/r49_1100_400_unitemis_2235g4.rtp'];
rtpwrite(rtp49x,hdr2,ha,pd,pa);
fin = rtp49x;
if exist('sym_ufo') delete('sym_ufo'); end
unix([' ln -s ' ufo  ' sym_ufo '])
command = [ SARTAEXE ' fin=' rtp49x ' fout=sym_ufo  > ' ...
            '/home/chepplew/logs/sarta/sar_stdout_unp.txt'];
unix(command)
[hdu hau pdu pau] = rtpread('sym_ufo'); 
sbc = rad2bt(hdr2.vchan, pdu.rcalc);
sbc_mn = nanmean(sbc,2);


kpath='/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2016_Mar2018/';

fnlst = dir([kpath 'RAD1013_unitemis/convolved_kcarta_RAD1013_*_radiances.mat']);

% reorder these in profile number order
  for i=1:49 fnparts{i} = strsplit(fnlst(i).name, '_'); 
    prfnums(i) = str2num(cell2mat(fnparts{i}(4)));
  end
  [B IB] = sort(prfnums);  
  krc = []; 
  for i=IB 
    x   = load([fnlst(i).folder '/' fnlst(i).name]);
    krc = [krc x.rcris_all];                 % [2235 x 392] 49 profs x 8 angles
  end

krc_nadir = krc(:,1:8:392);
kbc       = rad2bt(hdr2.vchan, krc_nadir);
kbc_mn    = nanmean(kbc,2);
bias_sd   = nanstd(kbc - sbc,0,2);


figure(1);clf;plot(hdr2.vchan, kbc_mn,'-', hdr2.vchan, sbc_mn,'-')
figure(1);clf;plot(hdr2.vchan, kbc_mn - sbc_mn,'-', hdr2.vchan, bias_sd,'-')
  axis([640 2600 -0.9 0.9]); grid on; legend('Mean Bias','Std.Dev of Bias');
  xlabel('Wavenumber cm^{-1}');ylabel('kCARTA minus SARTA (K)');
  title('kCARTA minus SARTA mean difference')
  
