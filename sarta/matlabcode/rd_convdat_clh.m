function [ctrans freq idchan] = rd_convdat_clh(filename)

% function [ctrans freq idchan] = rd_convdat_clh(filename)
% Read a conv trans data file
% Assumes data file is IEEE big endian

% Last updated: Scott Hannon, 20 July 2000, open big endian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fname=input('Enter name of file to read: ','s');
disp(['Reading convdata: ' filename]);
fname = filename;

% Open the input file
fid=fopen(fname,'r','ieee-be');

% Read the number of layers, gases, angles, and channels
ifm=fread(fid,1,'integer*4');
nang=fread(fid,1,'integer*4');
nlay=fread(fid,1,'integer*4');
ngas=fread(fid,1,'integer*4');
nchan=fread(fid,1,'integer*4');
ifm=fread(fid,1,'integer*4');


% Read the gas IDs
ifm=fread(fid,1,'integer*4');
gasid=fread(fid,ngas,'integer*4');
ifm=fread(fid,1,'integer*4');


% Read the angle secants
ifm=fread(fid,1,'integer*4');
angsec=fread(fid,nang,'real*4');
ifm=fread(fid,1,'integer*4');


% Read the comment
ifm=fread(fid,1,'integer*4');
junk=fread(fid,ifm,'char')';
comment=setstr(junk)
ifm=fread(fid,1,'integer*4');


% Read the temperature
ifm=fread(fid,1,'integer*4');
temp=fread(fid,nlay,'real*4');
ifm=fread(fid,1,'integer*4');


% Read the gas amounts
amount=zeros(nlay,ngas);
for igas=1:ngas
   ifm=fread(fid,1,'integer*4');
   amount(:,igas)=fread(fid,nlay,'real*4');
   ifm=fread(fid,1,'integer*4');
end


% Read the convolved transmittances
freq=zeros(1,nchan);
idchan=zeros(1,nchan);
res=zeros(1,nchan);
rnfwhm=zeros(1,nchan);
ctrans=zeros(nang,nlay,ngas,nchan);
for ichan=1:nchan
   ifm=fread(fid,1,'integer*4');
   freq(ichan)=fread(fid,1,'real*4');
   idchan(ichan)=fread(fid,1,'integer*4');
   res(ichan)=fread(fid,1,'real*4');
   rnfwhm(ichan)=fread(fid,1,'real*4');
   for igas=1:ngas
      for ilay=1:nlay
         ctrans(:,ilay,igas,ichan)=fread(fid,nang,'real*4'); 
      end
   end
   ifm=fread(fid,1,'integer*4');
end

% Close input file
ifm=fclose(fid);

% Clean up unwanted variables
clear ifm igas iang ichan
whos ctrans freq idchan res

%%% end of program %%%
