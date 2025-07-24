function wrt_concoef(outcoefname, idchan, fchan, coef);

% function wrt_concoef(outcoefname, idchan, fchan, coef);
%
% Write out continuum coefficients to fortran binary data file
%
% Input:
%    outcoefname : string, name of file to create
%    idchan      : integer(1,nchan) channel ID
%    fchan       : real(1,nchan) channel center frequency
%    coef        : real(ncoef,nlay,nchan) con coef
%

% Update: 05 Jan 2005, Scott Hannon - get array size from input not hardcoded
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check input
d = size(idchan);
if (length(d) ~= 2 | min(d) ~= 1)
   error('idchan must be a [1 x nchan] vector')
end
nchan = max(d)

d = size(fchan);
if (length(d) ~= 2 | min(d) ~= 1)
   error('fchan must be a [1 x nchan] vector')
end
ic = max(d);
if (ic ~= nchan)
   error('idchan and fchan must be the same length')
end

d = size(coef);
ncoef = d(1)
nlay = d(2)
ic = d(3);
if (ic ~= nchan)
   error('coef must have the same nchan as idchan and fchan')
end


% Record size for output file
irecsize=round( 4*( 2 + nlay*ncoef) );

fid=fopen(outcoefname,'w','ieee-be');
for ic=1:nchan

   ii=fwrite(fid,irecsize,'integer*4'); % start of record mark
   ii=fwrite(fid,idchan(ic),'integer*4'); % chan ID number
   ii=fwrite(fid,fchan(ic),'real*4'); % chan freq
   for il=1:nlay
      ii=fwrite(fid,coef(:,il,ic),'real*4');
   end
   ii=fwrite(fid,irecsize,'integer*4'); % end of record mark
end


% close output file
ii=fclose(fid);

%%% end of function %%%
