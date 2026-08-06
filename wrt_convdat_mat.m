function [iok]=wrt_convdat_mat(outname0, nang, nlay, ngas, nchan, gasid, ...
   secang, comment, temp, amount, freq, idchan, res, rnfwhm, ctrans, mfiles);

% function [iok]=wrt_convdat(outname0, nang, nlay, ngas, nchan, gasid, ...
%    secang, comment, temp, amount, freq, idchan, res, rnfwhm, ctrans);
%
% Function to write a FORTRAN format convolved layer-to-space transmittance
% binary data file for some profile.
%
% Input:
%    outname0 = string, name of file to create
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

outname = [outname0 '.mat'];

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
if exist(outname)
 outname
 %error('Error opening output file ... already exists')
 disp('warning output file ... already exists')
end

gasID = gasid(1:ngas);
secangs = secang(1:nang);
T = temp(1:nlay);
Q = amount(1:nlay,1:ngas);

L2Sdata.f  = freq(1:nchan);
L2Sdata.idchan = idchan(1:nchan);
L2Sdata.res    = res(1:nchan);
L2Sdata.rnfwhm = rnfwhm(1:nchan);
L2Sdata.ctrans = ctrans(:,1:nchan);

strall = [' nang nlay ngas nchan gasID secangs comment T Q L2Sdata mfiles'];

saver = ['save ' outname ' ' strall];
eval(saver);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iok=1;

lser = ['!ls -lt ' outname];
eval(lser)

iok=0;
thedir = dir(outname);
if length(thedir) == 1
  iok = 1;
  if thedir.bytes > 0
    iok = 2;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%% end of function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
