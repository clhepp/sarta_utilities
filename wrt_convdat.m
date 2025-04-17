function [iok]=wrt_convdat(outname, nang, nlay, ngas, nchan, gasid, ...
   secang, comment, temp, amount, freq, idchan, res, rnfwhm, ctrans);

% function [iok]=wrt_convdat(outname, nang, nlay, ngas, nchan, gasid, ...
%    secang, comment, temp, amount, freq, idchan, res, rnfwhm, ctrans);
%
% Function to write a FORTRAN format convolved layer-to-space transmittance
% binary data file for some profile.
%
% Input:
%    outname = string, name of file to create
%    nang = integer, number of angles
%    nlay = integer, number of layers
%    ngas = integer, number of gases
%    nchan = integer, number of channels
%    gasid = integer (1 x ngas), HITRAN gas ID numbers
%    secang = real (1 x nang), angle secants
%    comment = string (40 char), whatever comment to include in output file
%    temp = real (1 x nlay), temperature
%    amount = real (nlay x ngas), gas amounts
%    freq = real (1 x nchan), channel center freqs
%    idchan = integer (1 x nchan), channel ID numbers
%    res = real (1 x nchan), channel resolution
%    rnfwhm = real (1 x nchan), channel min(left,right) of FWHM into wing
%    ctrans = real ((nang*nlay*ngas) x nchan), convolved l-to-s transmittance
%
% Output:
%    iok = integer, 0 if this function detects a problem, otherwise 1
%
% Comment: This routine can not test for the order of the data in ctrans,
%    so it's up to you to make sure it's right.  This version forces
%    the output to be written using big endian IEEE standard fortran
%    format binary (ie with head & tail 4byte integer record markers).

% Created by Scott Hannon, 4 August 1999
% Last updated by Scott Hannon, 31 July 2000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize iok as bad
iok=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check array sizes are consistent
%
ijunk=length(gasid);
if (ijunk ~= ngas)
   error('Error, unexpected array size for gasid')
end
%
ijunk=length(secang);
if (ijunk ~= nang)
   error('Error, unexpected array size for secang')
end
%
ijunk=length(temp);
if (ijunk ~= nlay)
   error('Error, unexpected array size for temp')
end
%
[irow, ijunk]=size(amount);
if (irow ~= nlay | ijunk ~= ngas)
   error('Error, unexpected array size for amount')
end
%
ijunk=length(freq);
if (ijunk ~= nchan)
   error('Error, unexpected array size for freq')
end
%
ijunk=length(idchan);
if (ijunk ~= nchan)
   error('Error, unexpected array size for idchan')
end
%
ijunk=length(res);
if (ijunk ~= nchan)
   error('Error, unexpected array size for res')
end
%
ijunk=length(rnfwhm);
if (ijunk ~= nchan)
   error('Error, unexpected array size for rnfwhm')
end
%
[irow, ijunk]=size(ctrans);
if (irow ~= nang*nlay*ngas | ijunk ~= nchan)
   error('Error, unexpected array size for ctrans')
end


%%%%%%%%%%%%%%%%%%%%%%
% Open the output file
fid=fopen(outname,'w','ieee-be');
if (fid == -1)
   error('Error opening output file')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write header record: nang, nlay, ngas, nchan
ifm=4*4;
ic=fwrite(fid,ifm,'integer*4');
% Check first thing written to output file
if (ic ~= 1)
   error('Error writing to output file')
end
%
ic=fwrite(fid,nang,'integer*4');
ic=fwrite(fid,nlay,'integer*4');
ic=fwrite(fid,ngas,'integer*4');
ic=fwrite(fid,nchan,'integer*4');
%
%%%
%% old FORTRAN ordering
%ic=fwrite(fid,nlay,'integer*4');
%ic=fwrite(fid,ngas,'integer*4');
%ic=fwrite(fid,nang,'integer*4');
%ic=fwrite(fid,nchan,'integer*4');
%%%
%
ic=fwrite(fid,ifm,'integer*4');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write header record: gas IDs
ifm=ngas*4;
ic=fwrite(fid,ifm,'integer*4');
ic=fwrite(fid,gasid(1:ngas),'integer*4');
ic=fwrite(fid,ifm,'integer*4');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write header record: angle secants
ifm=nang*4;
ic=fwrite(fid,ifm,'integer*4');
ic=fwrite(fid,secang(1:nang),'real*4');
ic=fwrite(fid,ifm,'integer*4');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write header record: title
% Pad title/comment to make sure it is at least 40 characters
junk=[comment '                                        '];
ifm=40;
ic=fwrite(fid,ifm,'integer*4');
ic=fwrite(fid,junk(1:40),'char');
ic=fwrite(fid,ifm,'integer*4');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write header record: temperature
ifm=nlay*4;
ic=fwrite(fid,ifm,'integer*4');
ic=fwrite(fid,temp(1:nlay),'real*4');
ic=fwrite(fid,ifm,'integer*4');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write header records: gas amounts
% One record per gas
ifm=nlay*4;
for igas=1:ngas
   ic=fwrite(fid,ifm,'integer*4');
   ic=fwrite(fid,amount(1:nlay,igas),'real*4');
   ic=fwrite(fid,ifm,'integer*4');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write main data records: freq, idchan, res, rnfwhm, ctrans
% One record per channel
%
% Loop over the channels
ifm=(4 + irow)*4;
for ichan=1:nchan
   ic=fwrite(fid,ifm,'integer*4');
   ic=fwrite(fid,freq(ichan),'real*4');
   ic=fwrite(fid,idchan(ichan),'integer*4');
   ic=fwrite(fid,res(ichan),'real*4');
   ic=fwrite(fid,rnfwhm(ichan),'real*4');
   ic=fwrite(fid,ctrans(:,ichan),'real*4');
   ic=fwrite(fid,ifm,'integer*4');   
end
% Check last thing written to output file
if (ic ~= 1)
   error('Error writing to output file')
end

ijunk=fclose(fid);

iok=1;

%%%%%%%%%%%%%%%%%%%%%%%%% end of function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
