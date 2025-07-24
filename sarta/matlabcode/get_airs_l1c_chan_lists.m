% SARTA for AIRS L1C (2645)
%
% get_airs_l1c_chan_lists.m
% Routine to extend Scotts AIRS L1b channel lists to L1C 
% Load up Scott's channel lists for L1B sets 1 to 7.

addpath /asl/packages/airs_decon/source         % seq_match

a1 = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2378_list_set1');  % a1[1229x2]
a2 = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2378_list_set2');  % a1[303x2]
a3 = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2378_list_set3');  % a1[332x2]
a4 = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2378_list_set4');  % a1[67x2]
a5 = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2378_list_set5');  % a1[208x2]
a6 = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2378_list_set6');  % a1[158x2]
a7 = load('/home/chepplew/gitLib/ftc_dev/chanLists/airs_2378_list_set7');  % a1[81x2]
x=load('/home/chepplew/myLib/data/airs_f.mat');

% Check total number channels = 2378.
l1b_nchan = size(a1,1) + size(a2,1) + size(a3,1) + size(a4,1) + ...
            size(a5,1) + size(a6,1) + size(a7,1);
disp(['L1b channel total: ' num2str(l1b_nchan)]);


% The L1B set contains overlapping channels so need to apply unique operator
[xi xj]      = seq_match( sort(unique(a1(:,2))), x.fairs);    % length(x) = 1223.
cset(1).ich  = sort(unique(a1(:,1)));
cset(1).ich  = cset(1).ich(xi);
cset(1).frq  = x.fairs(xj);

[xi xj]     = seq_match( sort(unique(a2(:,2))), x.fairs);    % length = 303
cset(2).ich  = sort(unique(a2(:,1)));
cset(2).ich  = cset(2).ich(xi);
cset(2).frq  = x.fairs(xj);

[xi xj]     = seq_match( sort(unique(a3(:,2))), x.fairs);    % length = 320
cset(3).ich  = sort(unique(a3(:,1)));
cset(3).ich  = cset(3).ich(xi);
cset(3).frq  = x.fairs(xj);

[xi xj]     = seq_match( sort(unique(a4(:,2))), x.fairs);    % length = 67
cset(4).ich  = sort(unique(a4(:,1)));
cset(4).ich  = cset(4).ich(xi);
cset(4).frq  = x.fairs(xj);

[xi xj]     = seq_match( sort(unique(a5(:,2))), x.fairs);    % length = 175
cset(5).ich  = sort(unique(a5(:,1)));
cset(5).ich  = cset(5).ich(xi);
cset(5).frq  = x.fairs(xj);

[xi xj]     = seq_match( sort(unique(a6(:,2))), x.fairs);    % lenfth = 157
cset(6).ich  = sort(unique(a6(:,1)));
cset(6).ich  = cset(6).ich(xi);
cset(6).frq  = x.fairs(xj);

[xi xj]     = seq_match( sort(unique(a7(:,2))), x.fairs);    % length = 74
cset(7).ich  = sort(unique(a7(:,1)));
cset(7).ich  = cset(7).ich(xi);
cset(7).frq  = x.fairs(xj);

% check number channels included
lx_nchan = 0;
for i = 1:7 lx_nchan = lx_nchan + size(cset(i).frq,1); end     % 2329
disp(['Number of L1b channels included in L1C set: ' num2str(lx_nchan)]);


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
  


