function [ichan, fchan, coef, info] = rdcoef_wrt_netcdf(set, lmerged, fname);

% function [ichan, fchan, coef, info] = rdcoef_wrt_netcdf(set, lmerged, fname);
%
% Reads in binary FORTRAN coefficient data file and writes a netcdf file
%
% Input:
%    set     = {integer} set number [1:7} and 8=CO2/SO2/HNO3/NH3 4term, 
%              9=OPTRAN, 10=N2O 7term, 11=therm/nonLTE 12=CO2/SO2/HNO3 5term
%              13=HDO 11term
%    lmerged = {integer} 0=raw (not merged), 1=merged with water continuum
%    fname   = {string} name of binary FORTRAN data file to read
%
% Output:
%    ichan   = [nchan x 1] channel ID (JPL numbering)
%    fchan   = [nchan x 1] center frequency (wavenumber)
%    coef    = [nchan x nlay x ncoef] coeffcients
%    info    = {structure} various info
%
% Warning: if set=8 then info.gasid=2 even if SO2(9) or HNO3(12)
%

% 06 January 2004, Scott Hannon
% 02 Jan 2008, S.Hannon - add set and lcon to info; fix merged OPTRAN
%    reader; add set10 for N2O
% 09 Jan 2008, S.Hannon - add set for therm/nonLTE
% 15 May 2008. S.Hannon - add prompt for ncoef for nonlTE
% 29 Feb 2016. C. Hepplewhite - setX (X=1:8): nceof new values.
%        (use w/ /asl/data/sarta_database/Data_AIRS_apr08/Coef/setX_m130x370.dat)
% 10 May 2018 C Hepplewhite  Add NH3 set 8 4-term
%  1 Feb 2019 C Hepplewhite  Add set 13 for HDO 11term
%  1 Aug 2022 C Hepplewhite  Add capability to ready stacked xnte (14-coefficient).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify all allowed set numbers
allset  = [ 1  2  3  4  5  6  7  8  9 10 11 12 13];

if ( length(intersect(allset, set)) ~= 1 )
   allset
   set
   error('Unsupported set')
end

if (lmerged ~= 0 & lmerged ~= 1)
   lmerged
   error('lmerged must be zero (not merged) or 1 (merged)')
end

%command = ['stat -Lc %s ' fname];
%[a junk] = system(command);
%junk = str2num(junk);
junk = dir(fname);

if (length(junk) == 0)
   fname
   error('Did not find file')
end
nbytes = junk.bytes;
%disp(['nbytes: ' nbytes]);
fprintf(1,'nbytes : %i\n',nbytes)
clear junk


% Specify data for set1
if (set == 1)
   info.set=1;
   info.nlay=100;
   info.ncoef=24;                    % was 24=(8+11+5)
   info.nbreakout=3;
   info.gasid=[-2 1 3];
   info.ncoefeach=[8 11 5];
   info.startind=[1 9 20];
   if (lmerged == 1)
      info.lcon=1;
      ncon=7;
      info.ncoef=info.ncoef + ncon;
      info.startind=info.startind + ncon;
   else
      info.lcon=0;
   end
   xrecsize=4*(1 + 1 + info.nlay*info.ncoef);
   junk=nbytes/(xrecsize+8);
   info.nchan=round(junk);
   if (abs(junk - info.nchan) > 0.0001)
     error('unexpected file size') 
     info
     xrecsize
     junk
   end
end


% Specify data for set2
if (set == 2)
   info.set=2;
   info.nlay=100;
   info.ncoef=29;                 % was 29=(8+10+11)
   info.nbreakout=3;
   info.gasid=[-2 1 3];
   info.ncoefeach=[8 10 11];
   info.startind=[1 9 19];
   if (lmerged == 1)
      info.lcon=1;
      ncon=7;
      info.ncoef=info.ncoef + ncon;
      info.startind=info.startind + ncon;
   else
      info.lcon=0;
   end
   xrecsize=4*(1 + 1 + info.nlay*info.ncoef);
   junk=nbytes/(xrecsize+8);
   info.nchan=round(junk);
   if (abs(junk - info.nchan) > 0.0001)
     error('unexpected file size') 
     info
     xrecsize
     junk
   end
end


% Specify data for set3
if (set == 3)
   info.set=3;
   info.nlay=100;
   info.ncoef=28;                 % 28=(8+9+11) + 7 continuum;
   info.nbreakout=3;
   info.gasid=[-2 6 1];
   info.ncoefeach=[8 9 11];
   info.startind=[1 9 18];
   if (lmerged == 1)
      info.lcon=1;
      ncon=7;
      info.ncoef=info.ncoef + ncon;
      info.startind=info.startind + ncon;
   else
      info.lcon=0;
   end
   xrecsize=4*(1 + 1 + info.nlay*info.ncoef);
   junk=nbytes/(xrecsize+8);
   info.nchan=round(junk);
   if (abs(junk - info.nchan) > 0.0001)
     error('unexpected file size') 
     info
     xrecsize
     junk
   end
end


% Specify data for set4
if (set == 4)
   info.set=4;
   info.nlay=100;
   info.ncoef=38;              % 38=(11+11+3+13);
   info.nbreakout=4;
   info.gasid=[-2 5 3 1];
   info.ncoefeach=[11 11 3 13];
   info.startind=[1 12 23 26];
   if (lmerged == 1)
      info.lcon=1;
      ncon=7;
      info.ncoef=info.ncoef + ncon;
      info.startind=info.startind + ncon;
   else
      info.lcon=0;
   end
   xrecsize=4*(1 + 1 + info.nlay*info.ncoef);
   junk=nbytes/(xrecsize+8);
   info.nchan=round(junk);
   if (abs(junk - info.nchan) > 0.0001)
     error('unexpected file size') 
     info
     xrecsize
     junk
   end
end


% Specify data for set5
if (set == 5)
   info.set=5;
   info.nlay=100;
   info.ncoef=15;            % 15 = (11+3+1);
   info.nbreakout=3;
   info.gasid=[-2 1 3];
   info.ncoefeach=[11 3 1];
   info.startind=[1 12 15];
   if (lmerged == 1)
      info.lcon=1;
      ncon=7;
      info.ncoef=info.ncoef + ncon;
      info.startind=info.startind + ncon;
   else
      info.lcon=0;
   end
   xrecsize=4*(1 + 1 + info.nlay*info.ncoef);
   junk=nbytes/(xrecsize+8);
   info.nchan=round(junk);
   if (abs(junk - info.nchan) > 0.0001)
     error('unexpected file size') 
     info
     xrecsize
     junk
   end
end


% Specify data for set6
if (set == 6)
   info.set=6;
   info.nlay=100;
   info.ncoef=16;                  % 16=(8+7+1);
   info.nbreakout=3;
   info.gasid=[-2 1 3];
   info.ncoefeach=[8 7 1];
   info.startind=[1 9 16];
   if (lmerged == 1)
      info.lcon=1;
      ncon=7;
      info.ncoef=info.ncoef + ncon;
      info.startind=info.startind + ncon;
   else
      info.lcon=0;
   end
   xrecsize=4*(1 + 1 + info.nlay*info.ncoef);
   junk=nbytes/(xrecsize+8);
   info.nchan=round(junk);
   if (abs(junk - info.nchan) > 0.0001)
     error('unexpected file size') 
     info
     xrecsize
     junk
   end
end


% Specify data for set7
if (set == 7)
   info.set=7;
   info.nlay=100;
   info.ncoef=22;                 % 22=(8+13+1);
   info.nbreakout=3;
   info.gasid=[-2 1 3];
   info.ncoefeach=[8 13 1];
   info.startind=[1 9 22];
   if (lmerged == 1)
      info.lcon=1;
      ncon=7;
      info.ncoef=info.ncoef + ncon;
      info.startind=info.startind + ncon;
   else
      info.lcon=0;
   end
   xrecsize=4*(1 + 1 + info.nlay*info.ncoef);
   junk=nbytes/(xrecsize+8);
   info.nchan=round(junk);
   if (abs(junk - info.nchan) > 0.0001)
     error('unexpected file size') 
     info
     xrecsize
     junk
   end
end


% Specify data for set8 (perturbed CO2/SO2/HNO3)
if (set == 8)
   info.set=8;
   info.nlay=100;
   info.nbreakout=1;
   info.gasid=input('Enter gas ID (CO2=2, SO2=9, NH3=11, HNO3=12)');
%   info.gasid=[2]; % might be 2, 9, 11 or 12
   switch info.gasid
     case 2
       info.ncoef=4;              % was 4 for co2,so2
     case 9
       info.ncoef=4;
     case 11
       info.ncoef=4;
     case 12
       info.ncoef=4;
     otherwise
       disp('invalid gas id\n');
   end
   info.ncoefeach=[4];
   info.startind=[1];
   xrecsize=4*(1 + 1 + info.nlay*info.ncoef);
   junk=nbytes/(xrecsize+8);
   info.nchan=round(junk);
   if (abs(junk - info.nchan) > 0.0001)
      error('unexpected file size') 
   end
end


% Specify data for set9 (OPTRAN water)
if (set == 9)
   info.set=9;
   info.nlay=300;
   info.ncoef=9;
   info.nbreakout=1;
   info.gasid=[1];
   info.ncoefeach=[9];
   info.startind=[1];
   xrecsize=4*(1 + 1 + info.nlay*info.ncoef); % excludes record markers
   if (lmerged == 1)
      info.navgpred=4;
      info.avgpred=zeros(info.nlay,info.navgpred);
      hsize=4*(1 + info.navgpred)*(2 + info.nlay); % includes record markers
%wrong      junk=(nbytes-hsize)/xrecsize;
      junk=(nbytes-hsize)/(xrecsize+8);
   else
      junk=nbytes/(xrecsize+8);
   end
   info.nchan=round(junk);
   if (abs(junk - info.nchan) > 0.0001)
      error('unexpected file size') 
   end
end


% Specify data for set10 (perturbed N2O)
if (set == 10)
   info.set=10;
   info.nlay=100;
   info.ncoef=7;
   info.nbreakout=1;
   info.gasid=[4];
   info.ncoefeach=[4];
   info.startind=[1];
   xrecsize=4*(1 + 1 + info.nlay*info.ncoef);
   junk=nbytes/(xrecsize+8);
   info.nchan=round(junk);
   if (abs(junk - info.nchan) > 0.0001)
      error('unexpected file size') 
   end
end


% Specify data for set11 (therm/nonLTE)
if (set == 11)
   info.set=11;
   isnte=input('Is this non-LTE coefs? (0=no it is therm, 1=yes): ');
   if (isnte == 1)
      info.comment='daytime non-LTE coefs';
      junk = input('Enter number of nonLTE coefs (6 or 7 or 14) : ');
      info.ncoef = round(junk);
   else
      info.comment='reflected downwelling thermal F factor coefs';
      info.ncoef=6;
   end
   info.ncoefeach=info.ncoef;
   info.nlay=1;
   info.nbreakout=1;
   info.gasid=[-9999];
   info.startind=[1];
   xrecsize=4*(1 + 1 + info.nlay*info.ncoef);
   junk=nbytes/(xrecsize+8);
   info.nchan=round(junk);
   if (abs(junk - info.nchan) > 0.0001)
      error('unexpected file size') 
   end
end


% Specify data for set12 (5term CO2 OR 4term)
if (set == 12)
   info.set=12;
   info.nlay=100;
%   info.ncoef=4;
   info.ncoef=5;
   info.nbreakout=1;
info.gasid=input('Enter gas ID (CO2=2, SO2=9, HNO3=12)');
%   info.gasid=[2]; % might be 2, 9, or 12
   info.ncoefeach=[5];
   info.startind=[1];
   xrecsize=4*(1 + 1 + info.nlay*info.ncoef);
   junk=nbytes/(xrecsize+8);
   info.nchan=round(junk);
   if (abs(junk - info.nchan) > 0.0001)
      error('unexpected file size') 
   end
end

% Specify data for set 13 (11term HDO)
if ( set == 13)
  info.set = 13;
  info.nlay = 100;
  info.ncoef=11;
  info.nbreakout=1;
  info.gasid=162;         % AFGL code for HD^(16)O
  info.ncoefeach=[];
  info.startind=[];
  xrecsize=4*(1 + 1 + info.nlay*info.ncoef);
  junk=nbytes/(xrecsize+8);
  info.nchan=round(junk);
  if (abs(junk - info.nchan) > 0.0001)
     error('unexpected file size') 
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code below should not require modifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open output file
fid=fopen(fname,'r','ieee-be');


%%%
% Read merged OPTRAN header
if (set == 9 & lmerged == 1)
   [ifm,count]=fread(fid,1,'integer*4');
   if (ifm ~= 4*info.nlay)
      error('Unexpected number of elements in first record of optran header')
   end
   [info.ogrid,count]=fread(fid,info.nlay,'real*4');
   [ifm,count]=fread(fid,1,'integer*4');
   for ip=1:info.navgpred
      [ifm,count]=fread(fid,1,'integer*4');
      [info.avgpred(:,ip),count]=fread(fid,info.nlay,'real*4');
      [ifm,count]=fread(fid,1,'integer*4');
   end
end
%%%

% Dimension output arrays
ichan=zeros(info.nchan,1);
fchan=zeros(info.nchan,1);
coef=zeros(info.nchan,info.nlay,info.ncoef);

% Loop over the channels
for ic=1:info.nchan
   % Read FORTRAN start-of-record marker
   [ifm,count]=fread(fid,1,'integer*4');
   if (ifm ~= xrecsize)
      ifm
      xrecsize
      error('Unexpected coef data record size')
   end

   % Read data for this channel
   ichan(ic)=fread(fid,1,'integer*4');
   fchan(ic)=fread(fid,1,'real*4');
   coef(ic,:,:)=fread(fid,[info.ncoef,info.nlay],'real*4')'; %'

   % Read FORTRAN end-of-record marker
   ifm=fread(fid,1,'integer*4');
end % end of loop over bands

fclose(fid);
% =====================================================
% Write netcdf file
% ====================================================

junk = strsplit(fname,{'.'});
fnout = [junk{1} '.nc'];

nchan = info.nchan;
ncoef = info.ncoef;
nlay  = info.nlay;

nccreate(fnout,"ichan","Dimensions",{"r",nchan}); %,"Format","Classic")
nccreate(fnout,"fchan","Dimensions",{"r",nchan}); %,"Format","Classic")
nccreate(fnout,"coef","Dimensions",{"r",nchan,"l",nlay,"c",ncoef}); %,...
%   "Format","Classic")

ncwrite(fnout,"ichan",ichan);
ncwrite(fnout,"fchan",fchan);
ncwrite(fnout,"coef",coef);

% for OPTRAN need header values: ogrid [300x1], avgpred [300x4]
if(set == 9)
  nav   = info.navgpred;
  onlay = info.nlay;
  nccreate(fnout,"ogrid","Dimensions",{"n",onlay});
  nccreate(fnout,"avgpred","Dimensions",{"n",onlay,"o",nav});
  ncwrite(fnout,"ogrid", info.ogrid);
  ncwrite(fnout,"avgpred",info.avgpred);
end

%%% end of function %%%
