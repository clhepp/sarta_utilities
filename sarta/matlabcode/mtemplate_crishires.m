% CrIS hires version
% MATLAB script to convolve fast model layer-to-space transmittance data
% This version is for the Oct 2009 full OPD (hi-res) CrIS fast model
% production with four guard channels at each end of each band.  Note
% variable "ichan" is channel index but not channel ID.
%
% Basic outline of this script:
%    Run readkc3 to read in KCARTA data and create temporary output file
%    Read in the header of the temporary file
%    Figure out how many sets of 100-layer optical depths there are
%    Loop over the sets
%       Read in all layer-to-space optical depths for current set
%       Loop over angles
%          Convert nadir optical depth to transmittance for current angle
%          Convolve using sconv2
%       End loop on angles
%    End loop on sets
%    Save convolved trans
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The dummy X strings below should be replaced with the appropriate values
typeband='XTYPEBANDX';
outname='XMOUTNAMEX';
koutname='XKOUTNAMEX';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('Starting to process data ------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The variables below should be set by hand
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%
% nlay: Number of layers (per set of mixed paths)
nlay=100;


%%%%%%%
% secang: Angle secants (without sun)
secang=[1.000  1.190  1.410  1.680  1.990  2.370];


%%%%%%%
% secang_sun: Extra angle secants for shortwave
secang_sun=[2.840  3.470  4.300  5.420  6.940  9.020];


%%%%%%%
% tempfile: Name of the temporary file to be created by readkc3
tempfile='xxx.dat';

%%%%%%%
% alltypeband: Supported typeband
% allcfirst: First channel for each typeband
% allclast: Last channel for each typeband
alltypeband = { 'fowB1', 'fowB2', 'fmwB2', 'fcowB3'};
allband     = [       1,       2,       2,        3];
% hires CrIS with 0.8 cm OPD for all three bands
%%% without guard channels
%allcoffset  = [      72,      48,      48,       80];
%allcfirst   = [       1,     714,     714,     1579];
%allclast    = [     713,    1578,    1578,     2211];
%%% with four guard channels at each end
allcoffset  = [      68,      44,      44,       76];
allcfirst   = [       1,     722,     722,     1595];
allclast    = [     721,    1594,    1594,     2235];
%%%

%%% regular CrIS with 0.8, 0.4, 0.2 cm OPD
%allcoffset  = [      72,      24,      24,       20];
%allcfirst   = [       1,     714,     714,     1147];
%allclast    = [     713,    1146,    1146,     1305];
%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code below should not need modifying
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ii=strmatch(typeband,alltypeband);
if ( length(ii) > 1 | ii < 1)
   ii
   typeband
   alltypeband
   error(['invalid type_band: ' typeband]);
end


% Channel numbers
cfirst=allcfirst(ii);
clast=allclast(ii);
ichan=(cfirst:clast)'; %'
%iconv=ichan + allcoffset(ii);
iconv=(1:length(ichan))' + allcoffset(ii); %'
nchan=length(ichan);
clear alltypeband allcfirst allclast
clear cfirst clast
clear allcoffset

% Include sun angles
band = allband(ii);
if (band == 3)
   secang=[secang, secang_sun];
end
nang=length(secang);
clear secang_sun
clear allband

% Read in the kcarta output file and create a temporary output file
disp(' ')
disp(['Unchunking kcarta file ' koutname ' and creating temp file ' tempfile]);
readkc3(koutname,tempfile)
clear koutname

% Open and read in the header of the temporary file
disp(' ')
disp('Starting to processing data from temp file');
fid=fopen(tempfile,'r'); % Note: assumed written in native format
npts=fread(fid,1,'integer*4'); % number of freq points
ncol=fread(fid,1,'integer*4'); % number of paths
clear tempfile

% Read freq points
wnums=fread(fid,[npts,1],'real*8');

% Determine number of sets
nsets=ncol/nlay;
ijunk=mod(ncol,nlay);
if (ijunk ~= 0)
   ncol, nlay
   error('inappropriate value of ncol or nlay')
end
clear ijunk


% Add path for fftconv tools
addpath /asl/matlab/fftconv


% Convolve the layer-to-space transmittances
cris_ifpfile = ['cris12920B' int2str(band) 'hires'];
ctrans=zeros(nchan,ncol*nang);
indang=0:nang:round(nang*nlay - nang);  % indices for 1st angle of 1st set
% Loop over the sets
for iset=1:nsets
   disp(' ')
   disp(['Working on convolution set ' int2str(iset)]);
   % Read in all mono data for current set
   od=fread(fid,[npts,nlay],'real*4');
   %
   % Loop over angles
   for iang=1:nang
      disp(['doing angle number ' int2str(iang)])
      trans=exp(-od*secang(iang));
      %
      % Convolve
      [rout,fout] = xfconvkc(trans, cris_ifpfile, 'ham');
      ctrans(:,indang+iang)=rout(iconv,:);
   end % loop over angles
   %
   % Increment indang for next set
   indang=indang + nlay*nang;
   %
end % end loop over sets
fchan = fout(iconv)'; %'
clear iang iset indang npts ncol od trans rout fout wnums iconv

% Close temporary file
fclose(fid);
clear fid

% Save the results
disp(' ')
disp(['Saving convolved data to file ' outname]);
eval(['save ' outname ...
   ' ctrans fchan ichan secang nang nlay nsets nchan typeband'])

disp(' ')
disp('Finished processing data, quitting matlab')
exit

