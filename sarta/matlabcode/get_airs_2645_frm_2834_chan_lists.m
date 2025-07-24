function [iok] = get_airs_2645_frm_2378_chan_lists()

% Derive AIRS channel list sets for 2645 channels as used by AIRS L1C granules
% from Scott's 2834 lists
%
% Note that the frequencies of Scott's set are NOT EXACTLY the same for a given
%    channel than the standard set - so can't use intersect().

cd /home/chepplew/projects/sarta/airs/

addpath /asl/packages/airs_decon/source        % seq_match

% Get the nominal AIRS L1C (2645) channel set:
dpath      = '/asl/data/airs/L1C/2019/001/';
fn_l1c     = 'AIRS.2019.01.01.001.L1C.AIRS_Rad.v6.1.2.0.G19001103741.hdf';
junk       = hdfread([dpath fn_l1c], 'nominal_freq');
l1c.freq   = cell2mat(junk);
junk       = hdfread([dpath fn_l1c], 'ChanID');
l1c.chanID = cell2mat(junk);

fairs_2645 = Y1.fairs;
clear Y1;

% Get the 2834 airs channel set w/fake
junk = load('/asl/data/sarta_database/Data_AIRS_apr08/Coef/tunmlt_wcon_nte.txt');
airs_2834 = [junk(:,1), junk(:,2)];
clear junk;

% Get ordered indexes and matchups
[~, isx] = sort(l1c.freq);
[~, isz] = sort(airs_2834(:,2));

[iix, iiy] = seq_match(l1c.freq(isx), airs_2834(isz,2));


% Load the existing seven channel sets
airs_2834_set1 = load('/asl/data/sarta_database/Data_AIRS_apr08/List/list_set1');
airs_2834_set2 = load('/asl/data/sarta_database/Data_AIRS_apr08/List/list_set2');
airs_2834_set3 = load('/asl/data/sarta_database/Data_AIRS_apr08/List/list_set3');
airs_2834_set4 = load('/asl/data/sarta_database/Data_AIRS_apr08/List/list_set4');
airs_2834_set5 = load('/asl/data/sarta_database/Data_AIRS_apr08/List/list_set5');
airs_2834_set6 = load('/asl/data/sarta_database/Data_AIRS_apr08/List/list_set6');
airs_2834_set7 = load('/asl/data/sarta_database/Data_AIRS_apr08/List/list_set7');
airs_2834_optr = load('/asl/data/sarta_database/Data_AIRS_apr08/List/list_optran');


% Match up these with the standard AIRS set then assign ichan and fchan values.
%%% set1
[IX IY] = seq_match(sort(airs_2834_set1(ind1,2)),sort(fairs_2378));
if(length(IY) < length(ind1))
   missing = setdiff(ind1, IX);
   IY = sort([IY' IY(missing)]);
end
airs_2378_set1  = [IY', fairs_2378(IY)];
%%% set2
[IX IY] = seq_match(sort(airs_2834_set2(ind2,2)),sort(fairs_2378));
if(length(IY) < length(ind2))
   missing = setdiff(ind1, IX);
   IY = sort([IY' IY(missing)]);
end
airs_2378_set2  = [IY, fairs_2378(IY)];
%%% set3
[IX IY] = seq_match(sort(airs_2834_set3(ind3,2)),sort(fairs_2378));
if(length(IY) < length(ind3))
   missing = setdiff(ind3, IX);
   IY = sort([IY' IY(missing)]);
end
airs_2378_set3  = [IY, fairs_2378(IY)];
%%% set4
[IX IY] = seq_match(sort(airs_2834_set4(ind4,2)),sort(fairs_2378));
if(length(IY) < length(ind4))
   missing = setdiff(ind4, IX);
   IY = sort([IY' IY(missing)]);
end
airs_2378_set4  = [IY, fairs_2378(IY)];
%%% set5
[IX IY] = seq_match(sort(airs_2834_set5(ind5,2)),sort(fairs_2378));
if(length(IY) < length(ind5))
   missing = setdiff(ind5, IX);
   IY = sort([IY' IY(missing)']);
end
airs_2378_set5  = [IY', fairs_2378(IY)];
%%% set6
[IX IY] = seq_match(sort(airs_2834_set6(ind6,2)),sort(fairs_2378));
if(length(IY) < length(ind6))
   missing = setdiff(ind6, IX);
   IY = sort([IY' IY(missing)']);
end
airs_2378_set6  = [IY', fairs_2378(IY)];
%%% set7
[IX IY] = seq_match(sort(airs_2834_set7(ind7,2)),sort(fairs_2378));
if(length(IY) < length(ind7))
   missing = setdiff(ind7, IX);
   IY = sort([IY' IY(missing)']);
end
airs_2378_set7  = [IY', fairs_2378(IY)];
%%% Optran
[IX IY] = seq_match(sort(airs_2834_optr(indoptr,2)),sort(fairs_2378));
if(length(IY) < length(indoptr))
   missing = setdiff(indoptr, IX);
   IY = sort([IY' IY(missing)']);
end
airs_2378_optr  = [IY', fairs_2378(IY)];


%{
% Sanity check

figure(1);clf; 
  plot([1:2378],fairs_2378,'.',airs_2378_optr(:,1),airs_2378_optr(:,2),'o'); 
figure(2);clf;
   plot(airs_2834(:,1),airs_2834(:,2),'.',airs_2834_optr(:,1),airs_2834_optr(:,2),'o');
%}

% Write text files out
FH = fopen('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2378_list_set1','w');
for i =1:length(airs_2378_set1)
  fprintf(FH, '%4i\t%9.3f\n', airs_2378_set1(i,1), airs_2378_set1(i,2));
end
fclose(FH);
FH = fopen('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2378_list_set2','w');
for i =1:length(airs_2378_set2)
  fprintf(FH, '%4i\t%9.3f\n', airs_2378_set2(i,1), airs_2378_set2(i,2));
end
fclose(FH);
FH = fopen('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2378_list_set3','w');
for i =1:length(airs_2378_set3)
  fprintf(FH, '%4i\t%9.3f\n', airs_2378_set3(i,1), airs_2378_set3(i,2));
end
fclose(FH);
FH = fopen('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2378_list_set4','w');
for i =1:length(airs_2378_set4)
  fprintf(FH, '%4i\t%9.3f\n', airs_2378_set4(i,1), airs_2378_set4(i,2));
end
fclose(FH);
FH = fopen('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2378_list_set5','w');
for i =1:length(airs_2378_set5)
  fprintf(FH, '%4i\t%9.3f\n', airs_2378_set5(i,1), airs_2378_set5(i,2));
end
fclose(FH);
FH = fopen('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2378_list_set6','w');
for i =1:length(airs_2378_set6)
  fprintf(FH, '%4i\t%9.3f\n', airs_2378_set6(i,1), airs_2378_set6(i,2));
end
fclose(FH);
FH = fopen('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2378_list_set7','w');
for i =1:length(airs_2378_set7)
  fprintf(FH, '%4i\t%9.3f\n', airs_2378_set7(i,1), airs_2378_set7(i,2));
end
fclose(FH);
FH = fopen('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2378_list_optran','w');
for i =1:length(airs_2378_optr)
  fprintf(FH, '%4i\t%9.3f\n', airs_2378_optr(i,1), airs_2378_optr(i,2));
end
fclose(FH);
  
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Get AIRS solar flux spectrum
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FH = fopen('/asl/data/sarta_database/Data_AIRS_apr08/Solar/solar_m140x.txt','r');
iRow = 1;
while (~feof(FH)) 
    solData(iRow,:) = textscan(FH,'%f %f %f\n','CommentStyle','!');
    iRow = iRow + 1;
end
fclose(FH);

% Interpolate onto the airs_2378 grid
sol_airs_2378 = interp1(solData{2}, solData{3}, fairs_2378);

% write new solar data file:
FH = fopen('/asl/s1/chepplew/data/sarta_database/Data_AIRS_2378_400ppm/Solar/sol.txt','w');
for i = 1:length(sol_airs_2378)
  fprintf(FH,'%5i  %9.3f  %9.4f\n',i,fairs_2378(i),sol_airs_2378(i));
end
fclose(FH);

%    %%%%%%%%%%%%%%%%%%%%%%%%%%
%    tunmult_2378.txt all ones
%    %%%%%%%%%%%%%%%%%%%%%%%%%%
myOnes = ones(2378,7);
FH = fopen('/asl/s1/chepplew/data/sarta_database/Data_AIRS_2378_400ppm/Coef/tunmlt_2378.txt','w');
for i=1:2378
  fprintf(FH,'%5i  %9.3f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n', ...
          i,fairs_2378(i),myOnes(i,:));
end
fclose(FH);

   
