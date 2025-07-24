% SARTA for CHIRP (1702 channels)
%
% get_chirp_chan_lists.m
% Routine to match Scotts AIRS L1b channel lists to CHIRP grid 
% Load up Scott's ORIGINAL channel lists for AIRS L1B with fake.

cd /home/chepplew/projects/sarta

addpath /asl/matlib/h4tools                     % rtpread
addpath /asl/packages/airs_decon/source         % seq_match
addpath /asl/matlib/aslutil                     % rad2bt

% Get CHIRP spectrum
ccst=load('/asl/cris/ccast/sdr45_npp_MR/2018/250/CrIS_SDR_npp_s45_d20180907_t1130080_g116_v20a.mat');
frq_msr  = [ccst.vLW' ccst.vMW' ccst.vSW'];
rad_msr  = [ccst.rLW; ccst.rMW; ccst.rSW];
junk     = reshape(rad_msr,[], 9*30*45);
rad_msr_ham = hamm_app(double(junk));
ccst_bt  = real(rad2bt(frq_msr, rad_msr_ham));
ccst_btm = mean(ccst_bt,2);

% Get AIRS spectrum
SARTAEXE='/asl/packages/sartaV108/BinV201/sarta_apr08_m140_wcon_nte';
dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
             'REGR49_400ppm_H2016_Feb2020_AIRS2834_CHIRP/'];
rtpfile   = [dpath 'regr49_1100_400ppm_unitemiss.op.rtp'];
FIN = 'sarta_in.rtp';
FOUT = 'junk_out.rtp';
eval(['! ' SARTAEXE ' fin=' FIN ' fout=' FOUT ' > /home/chepplew/logs/sarta/sar_out.log']);

[hd ha pd pa] = rtpread('junk_out.rtp');
bbt = rad2bt(hd.vchan, mean(pd.rcalc,2));

% Load AIRS channel Lists, assign array of names (in order)
names = {'set1','set2','set3','set4','set5','set6','set7',...
         'hno3','n2o','co2_5term','so2','nh3','optran'};

a1  = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2834_list_set1');  % a1[1461x2]
a2  = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2834_list_set2');  % a1[325x2]
a3  = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2834_list_set3');  % a1[396x2]
a4  = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2834_list_set4');  % a1[85x2]
a5  = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2834_list_set5');  % a1[210x2]
a6  = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2834_list_set6');  % a1[217x2]
a7  = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2834_list_set7');  % a1[140x2]
a8  = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2834_list_hno3');  % [# x 2]
a9  = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2834_list_n2o');  % [# x 2]
a10 = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2834_list_co2_5term');  % [# x 2]
a11 = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2834_list_so2');  % [# x 2]
a12 = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2834_list_nh3');  % [# x 2]
a13 = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2834_list_optran');  % [754 x 2]

fairs = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2834_list_all');
x=load('/home/chepplew/myLib/data/chirp_1702_freq.mat');

% Check total number channels for lists 1..7 = 2378.
airs_nchan = size(a1,1) + size(a2,1) + size(a3,1) + size(a4,1) + ...
             size(a5,1) + size(a6,1) + size(a7,1);
disp(['AIRS channel total: ' num2str(airs_nchan)]);

% --------------- Sort lists 1..7 ---------------------
% The L1B set contains overlapping channels so need to apply unique operator
[xi xj]      = seq_match(unique(sort(a1(:,2))), x.vchan);    % length = 788
junk         = sort(unique(a1(:,1)));
cset(1).ich  = junk(xi,1);
cset(1).frq  = x.vchan(xj);

[xi xj]     = seq_match(unique(sort(a2(:,2))), x.vchan);    % length = 194
%cset(2).ich  = sort(unique(a2(:,1)));
cset(2).ich  = xj; % a2(xi,1);
cset(2).frq  = x.vchan(xj);

[xi xj]     = seq_match(unique(sort(a3(:,2))), x.vchan);    % length = 266
%cset(3).ich  = sort(unique(a3(:,1)));
cset(3).ich  = xj; % a3(xi,1);
cset(3).frq  = x.vchan(xj);

[xi xj]     = seq_match(unique(sort(a4(:,2))), x.vchan);    % length = 62
%cset(4).ich  = sort(unique(a4(:,1)));
cset(4).ich  = xj; % a4(xi,1);
cset(4).frq  = x.vchan(xj);

[xi xj]     = seq_match(unique(sort(a5(:,2))), x.vchan);    % length = 124
%cset(5).ich  = sort(unique(a5(:,1)));
cset(5).ich  = xj; % unique(a5(xi,1));
cset(5).frq  = x.vchan(xj);

[xi xj]     = seq_match(unique(sort(a6(:,2))), x.vchan);    % length = 125
%cset(6).ich  = sort(unique(a6(:,1)));
cset(6).ich  = xj; % unique(a6(xi,1));
cset(6).frq  = x.vchan(xj);

[xi xj]     = seq_match(unique(sort(a7(:,2))), x.vchan);    % length = 3
%cset(7).ich  = sort(unique(a7(:,1)));
cset(7).ich  = xj; % a7(xi,1);
cset(7).frq  = x.vchan(xj);

% check number channels included
lx_nchan = 0;
for i = 1:7 lx_nchan = lx_nchan + size(cset(i).frq,1); end     % 2329
disp(['Number of L1b channels included in L1C set: ' num2str(lx_nchan)]);

% Concatenate
vchan_all = [cset(1).frq' cset(2).frq' cset(3).frq' cset(4).frq' ...
             cset(5).frq' cset(6).frq' cset(7).frq' ];

% Find the gaps with the fake channels
f_all      = [cset(1).frq' cset(2).frq' cset(3).frq' cset(4).frq' ...
              cset(5).frq' cset(6).frq' cset(7).frq'];
[zi zj]    = setdiff(x.fairs, f_all);
gap(1).frq  = zi(1:21);
gap(2).frq  = zi(22:42);
gap(3).frq  = zi(43:63);
gap(4).frq  = zi(64:84);
gap(5).frq  = zi(85:237);
gap(6).frq  = zi(238:258);
gap(7).frq  = zi(259:287);
gap(8).frq  = zi(288:308);
gap(9).frq  = zi(309:331);      % 

% Now fill in the gaps - choose which set to apply the extra channels.
% gaps 1, 2 & 3 to set 1
% gaps 4 & 5    to set 2
% gaps 6, 7 & 8 to set 3
% gap  9        to set 6
[cset(1).fill sj] = sort([cset(1).frq' gap(1).frq' gap(2).frq' gap(3).frq']);         
                   % length() = 1286
[cset(2).fill sj] = sort([cset(2).frq' gap(4).frq' gap(5).frq']);         
                   % length() = 477
[cset(3).fill sj] = sort([cset(3).frq' gap(6).frq' gap(7).frq' gap(8).frq']);         
                   % length() = 401
[cset(4).fill sj] = sort([cset(4).frq']);
[cset(5).fill sj] = sort([cset(5).frq']);
[cset(6).fill sj] = sort([cset(6).frq' gap(9).frq']);
[cset(7).fill sj] = sort([cset(7).frq']);

[~, cset(1).cch]  = seq_match(cset(1).fill, x.fairs);   % no repeats
[~, cset(2).cch]  = seq_match(cset(2).fill, x.fairs);   % no repeats
[~, cset(3).cch]  = seq_match(cset(3).fill, x.fairs);   % no repeats
[~, cset(4).cch]  = seq_match(cset(4).fill, x.fairs);   % no repeats
[~, cset(5).cch]  = seq_match(cset(5).fill, x.fairs);   % no repeats
[~, cset(6).cch]  = seq_match(cset(6).fill, x.fairs);   % no repeats
[~, cset(7).cch]  = seq_match(cset(7).fill, x.fairs);   % no repeats

l1c_nchan = 0; l1c_ichan = [];
for i=1:7 
  l1c_nchan = l1c_nchan + size(cset(i).fill,2); 
  l1c_ichan = [l1c_ichan cset(i).cch'];
end
l1c_fall = [];
for i=1:7 l1c_fall = [l1c_fall cset(i).fill]; end
disp(['Number of channels in filled L1C set: ' num2str(l1c_nchan)]);

% I have an overlaps between sets 1&2 (2), sets 5&6 (11) and sets 6&7 (2).
junk = diff(sort(l1c_ichan));
numel(find(junk == 0))               % 15

[ii iy]   = intersect(cset(1).cch, cset(2).cch);     % 2
  cset(1).fill(iy) = [];
  cset(1).cch(iy)  = [];
[ii iy]   = intersect(cset(5).cch, cset(6).cch);     % 10
  cset(5).fill(iy) = [];
  cset(5).cch(iy)  = [];
[ii iy]   = intersect(cset(5).cch, cset(7).cch);     % 2
  cset(5).fill(iy) = [];
  cset(5).cch(iy)  = [];
[ii iy]   = intersect(cset(6).cch, cset(7).cch);     % 1
  cset(6).fill(iy) = [];
  cset(6).cch(iy)  = [];

l1c_nchan = 0; l1c_ichan = [];
for i=1:7 
  l1c_nchan = l1c_nchan + size(cset(i).fill,2); 
  l1c_ichan = [l1c_ichan cset(i).cch'];
end
l1c_fall = [];
for i=1:7 l1c_fall = [l1c_fall cset(i).fill]; end
disp(['Number of channels in final L1C set: ' num2str(l1c_nchan)]);  % 2645

% save data for use with ftc_dev fitting routines
savd = '/home/chepplew/gitLib/ftc_dev/chanLists/';
for i=1:7
  savf = sprintf('airs_2645_list_set%d',i);
  FH = fopen([savd savf],'w');
    fprintf(FH,'  %i\t %f\n', [cset(i).cch'; cset(i).fill]);
  fclose(FH);
end

% for the water continuum need full set in ascii text file)
savf = 'airs_2645_list_all';
FH = fopen([savd savf], 'w');
junk = [sort(l1c_ichan); sort(l1c_fall)];
fprintf(FH, '  %i\t %f\n', junk);
fclose(FH);
 

%{
fh1=figure(1);set(fh1,'resize','off');set(fh1,'Position',fh1.Position+[0 0 240 0]);
  plot(cset(1).ich, cset(1).frq,'.'); hold on;
  plot(cset(2).ich, cset(2).frq,'.');
  plot(cset(3).ich, cset(3).frq,'.');
  plot(cset(4).ich, cset(4).frq,'.');
  plot(cset(5).ich, cset(5).frq,'.');
  plot(cset(6).ich, cset(6).frq,'.');
  plot(cset(7).ich, cset(7).frq,'.');
  legend('1','2','3','4','5','6','7')


fh2=figure(2);set(gcf,'resize','off');set(fh2,'Position',fh2.Position+[0 0 240 0]);
figure(2);clf; plot(cset(1).cch, cset(1).fill ,'.');
  hold on;plot(cset(2).cch, cset(2).fill,'o');
  plot(cset(3).cch, cset(3).fill,'d');
  plot(cset(4).cch, cset(4).fill,'.');
  plot(cset(5).cch, cset(5).fill,'o');
  plot(cset(6).cch, cset(6).fill,'d');
  plot(cset(7).cch, cset(7).fill,'.');
  
figure(3);clf;plot([1:2645],x.fairs,'.', l1c_ichan, l1c_fall,'.')

%}
  
%{

% 2. Manually select channel lists by visual comparison with Scott's original

[ix{1} iy{1}]  = seq_match(sort(hd.vchan), sort(a1(:,2)));
[ix{2} iy{2}]  = seq_match(sort(hd.vchan), sort(a2(:,2)));
[ix{3} iy{3}]  = seq_match(sort(hd.vchan), sort(a3(:,2)));
[ix{4} iy{4}]  = seq_match(sort(hd.vchan), sort(a4(:,2)));
[ix{5} iy{5}]  = seq_match(sort(hd.vchan), sort(a5(:,2)));
[ix{6} iy{6}]  = seq_match(sort(hd.vchan), sort(a6(:,2)));
[ix{7} iy{7}]  = seq_match(sort(hd.vchan), sort(a7(:,2)));
figure(1);clf;plot(hd.vchan, bbt,'o-')
  hold on; 
  plot(hd.vchan(ix{1}), bbt(ix{1}),'r.')
  plot(hd.vchan(ix{2}), bbt(ix{2}),'g.')
  plot(hd.vchan(ix{3}), bbt(ix{3}),'r.')
  plot(hd.vchan(ix{4}), bbt(ix{4}),'g.')
  plot(hd.vchan(ix{5}), bbt(ix{5}),'r.')
  plot(hd.vchan(ix{6}), bbt(ix{6}),'g.')
  plot(hd.vchan(ix{7}), bbt(ix{7}),'m.')
xlim([650 700])

cset(1).ich = [1:543 681:702 709:710];
cset(1).fch = frq_msr(cset(1).ich);
%
cset(2).ich = [544:680 703:708 711:717];
cset(2).fch = frq_msr(cset(2).ich);
%
cset(3).ich = [718:1370];
cset(3).fch = frq_msr(cset(3).ich);
%
cset(4).ich = [1371:1443];
cset(4).fch = frq_msr(cset(4).ich);
%
cset(5).ich = [1444:1563];
cset(5).fch = frq_msr(cset(5).ich);
%
cset(6).ich = [1564:1691]; % [1564:1702];
cset(6).fch = frq_msr(cset(6).ich);
%
cset(7).ich = [];
cset(7).fch = frq_msr(cset(7).ich);

[jx{1} jy{1}]  = seq_match(frq_msr, cset(1).fch);
[jx{2} jy{2}]  = seq_match(frq_msr, cset(2).fch);
[jx{3} jy{3}]  = seq_match(frq_msr, cset(3).fch);
[jx{4} jy{4}]  = seq_match(frq_msr, cset(4).fch);
[jx{5} jy{5}]  = seq_match(frq_msr, cset(5).fch);
[jx{6} jy{6}]  = seq_match(frq_msr, cset(6).fch);
[jx{7} jy{7}]  = seq_match(frq_msr, cset(7).fch);
hold on;
  plot(frq_msr, ccst_btm+10, 'ko-')
  plot(frq_msr(jx{1}), ccst_btm(jx{1})+10,'r.')
  plot(frq_msr(jx{2}), ccst_btm(jx{2})+10,'g.')
  plot(frq_msr(jx{3}), ccst_btm(jx{3})+10,'r.')
  plot(frq_msr(jx{4}), ccst_btm(jx{4})+10,'g.')
  plot(frq_msr(jx{5}), ccst_btm(jx{5})+10,'r.')
  plot(frq_msr(jx{6}), ccst_btm(jx{6})+10,'g.')
  plot(frq_msr(jx{7}), ccst_btm(jx{7})+10,'c.')
%}

%{
% Write Channel Lists
outd = '/home/chepplew/gitLib/ftc_dev/chanLists/';
for i=1:7
  fn(i).name = ['chirp_1702_list_set' num2str(i)];
end
for i=1:7
FH = fopen([outd fn(i).name],'w');
  for k=1:length(cset(i).ich)
    fprintf(FH, '  %7i   %8.3f\n', cset(i).ich(k), cset(i).fch(k));
  end
fclose(FH);
%}

% =================================================
% Other breakouts

%
[junk, iix] = unique(a10(:,2));
aX = a10; clear a10; a10 = [ aX(iix,1), junk];

[ix{8} iy{8}]    = seq_match(sort(hd.vchan), sort(a8(:,2)));          % HNO3
[ix{9} iy{9}]    = seq_match(sort(hd.vchan), sort(a9(:,2)));          % N2O
[ix{10} iy{10}]  = seq_match(sort(hd.vchan), sort(a10(:,2)));         % CO2_5
[ix{11} iy{11}]  = seq_match(sort(hd.vchan), sort(a11(:,2)));         % SO2
[ix{12} iy{12}]  = seq_match(sort(hd.vchan), sort(a12(:,2)));         % NH3
[ix{13} iy{13}]  = seq_match(sort(hd.vchan), sort(a13(:,2)));         % Optran 649
%
figure(1);clf;plot(hd.vchan, bbt,'o-')
  hold on; 
  plot(hd.vchan(ix{9}), bbt(ix{9}),'r.')
hold on;
  plot(frq_msr, ccst_btm-10, 'ko-')

% Read off the CHIRP channels corresponding to AIRS for selected gas
%
cset(8).ich = [156:236 304:446 789:901 1273:1370];
cset(8).fch = frq_msr(cset(8).ich);
%
cset(9).ich = [743:867 1371:1458 1586:1646 1666:1691];
cset(9).fch = frq_msr(cset(9).ich);
%
cset(10).ich = [1:262 460:540 596:717 1419:1642];
cset(10).fch = frq_msr(cset(10).ich);
%
cset(11).ich = [663:763 843:951 1613:1672];
cset(11).fch = frq_msr(cset(11).ich);
%
cset(12).ich = [164:717 718:751 1009:1370];
cset(12).fch = frq_msr(cset(12).ich);
%
cset(13).ich = [170 193:195 204:206 210 217:219 229:231 236:237 240:242 ...
                247:250 255:257 265:267 283:284 286:288 306:308 322:323 ...
		326:328 330:331 347:348 356:358 368:370 376:378 382:383 ...
		416:418 437:440 442:444 479:481 485 487:489 524:525 ...
		682:683 693:694 708:709 ...
		721:724 730:731 737:739 743:744 751:753 760:763 765 ...
		772:773 779:781 787  790:791 794 806:809 812:815 824:826 ...
		842:845 848:852 860:861 866:868 872:879 884 887 ...
		894:895 898 901:905 910:911 916:919 922:924 932:934 937 941:943 947 ...
		948:1370];
cset(13).fch = frq_msr(cset(13).ich);

[jx{8} jy{8}]    = seq_match(frq_msr, cset(8).fch);
[jx{9} jy{9}]    = seq_match(frq_msr, cset(9).fch);
[jx{10} jy{10}]  = seq_match(frq_msr, cset(10).fch);
[jx{11} jy{11}]  = seq_match(frq_msr, cset(11).fch);
[jx{12} jy{12}]  = seq_match(frq_msr, cset(12).fch);
[jx{13} jy{13}]  = seq_match(frq_msr, cset(13).fch);

% Visual inspect plot
hold on
  plot(frq_msr(jx{9}), ccst_btm(jx{9})-10,'r.')

% Write channel list files
outd = '/home/chepplew/gitLib/ftc_dev/chanLists/';

% co2_5term: 10. optran: 13
nk = 9;
fn(nk).name = ['chirp_1702_list_' str2mat(names{nk})];
FH = fopen([outd fn(nk).name],'w');
  for k=1:length(cset(nk).ich)
    fprintf(FH, '  %7i   %8.3f\n', cset(nk).ich(k), cset(nk).fch(k));
  end
fclose(FH);
   

