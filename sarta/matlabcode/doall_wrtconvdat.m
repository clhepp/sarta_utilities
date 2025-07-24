function [ok]=doall_wrtconvdat(comment, pnums, rtpfile, matfile, outfile, lfcow2fow);

% function [ok]=doall_wrtconvdat(comment, pnums, rtpfile, matfile, outfile,lfcow2fow);
%
% Routine to merge the matlab convolved l-to-s trans data (created
% by the "runtrans" kcarta & matlab script) and the fortran profile
% data (from an RTP file) and write out a fortran binary data
% file with all the data in the format needed for use with our
% fast model regression programs.  It loops over all the regression
% profiles and writes a sperate output file for each profile.
%
% Input:
%    comment = {string (max 35 char)}, comment string which will be
%       included in the comment/title portion of the output files.
%       The profile number will be appended to the end as " p<#>".
%       Example: ['6 Aug 99, SRF 149 Aug99, fow' ' p13']
%    pnums = [1 x nprof] list of integer profile numbers {eg [1:48]}
%    rtpfile = {string} name of RTP file with regression data profile
%    matfile = {string} name prefix of matlab convolved trans files.
%       The file names are assumed to be:
%       <whatever><pnum>'.mat'
%    outfile = string, name prefix of fortran output file to create.
%       The output file names will be:
%       <outfile><pnum>'.dat'
%    lfcow2fow = OPTIONAL [1 x 1] convert fcowB3 input to fowB3
%       output? (0=no=default, 1=yes). When "yes", comment and outfile
%       should be for fowB3.
%
% Output:
%    ok = integer, 0 if the function detects a problem, otherwise 1
%
% Note: uses function "wrt_convdat.m".
%

% Original created by Scott Hannon, 6 August 1999
% Extensively revised to read in the profile data from RTP; April 2002
% Last updated: 9 May 2002, Scott Hannon - fix bug & check comment
% 13 Apr 2009, S.Hannon - update alltypeband for CrIS; add lfcow2fow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Info for supported typeband
alltypeband={'fowB1', 'fowB2', 'fmwB2', 'fcowB3'};
allngas=[       4,       4,       3,       5    ]; 
allgasids=zeros( length(alltypeband), max(allngas) );
allgasids(1,1:4)=[2, 3, 1, -2];     % fowB1
allgasids(2,1:4)=[2, 3, 1, -2];     % fowB2
allgasids(3,1:3)=[2, 6, 1];         % fmwB2
allgasids(4,1:5)=[2, 5, 3, 1, -2];  % fcowB3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ok=0;

if (nargin < 5 | nargin > 6)
   error('Unexpected number of input arguements')
end
if (nargin == 5)
   lfcow2fow = 0;
end

% Check comment
if (length(comment) > 35)
   error('comment string is too long; max allowed length is 35 char')
end


% Read the RTP profile data
addpath /asl/matlab/h4tools
[head, hattr, prof, pattr] = rtpread(rtpfile);
if (head.ptype ~= 1)
   error('rtpfile is not a "layers" profile')
end

% Loop over the profiles
for ip=1:length(pnums)
   ipnum=pnums(ip);
   disp(['processing profile ' int2str(ipnum)])

   %%%%%%%%%%%%%%%%%%%%%%%%%
   % Read in the matlab data
   % The matfile must contain:
   %       ctrans(nchan x [nang,nlay,nsets]) {layers bot to top}
   %       secang(nang)
   %       fchan(nchan)
   %       ichan(nchan)
   %       typeband
   %       nchan, nang, nlay, nsets
   mfile=[matfile int2str(ipnum) '.mat'];
   djunk=dir(mfile);
   [irow,icol]=size(djunk);
   if (irow ~= 1)
      error(['Error openning matlab file ' mfile])
   end
   clear djunk irow icol
   eval(['load ' mfile]);


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Check typeband and determine gas IDs
   ii=strmatch(typeband,alltypeband);
   if ( length(ii) > 1 | ii < 1)
      ii
      typeband
      alltypeband
      ipnum
      error(['invalid type_band: ' typeband]);
   end
if (ii ~= 4)
lfcow2fow = 0;
end
   ngas=allngas(ii);
   if (ngas ~= nsets)
     ngas
     nsets
     ipnum
     error('mismatch in expected ngas and matfile nsets')
   end
   gasids=allgasids(ii,1:ngas);


   %%%%%%%%%%%%%%%%%%
   % Load temperature
   ii=prof.nlevs(ipnum) - 1; % number of layers
   if (ii ~= nlay)
      ii
      nlay
      ipnum
      error('mismatch in number of layers')
   end
   temp=prof.ptemp(1:nlay,ipnum);


   %%%%%%%%%%%%%%%%%%
   % Load gas amounts
   amount=zeros(nlay,ngas);
   for ig=1:ngas

      ii=intersect( head.glist, abs(gasids(ig)) );
      if (length(ii) ~= 1)
         ipnum
         gasids(ig)
         head.glist
         error('required gas not found in RTP file')
      end
      ii=abs(gasids(ig));
      eval(['amount(:,ig)=prof.gas_' int2str(ii) '(1:nlay,ipnum);'])
   end
   % Convert amount from molecules/cm^2 to kmoles/cm^2
   amount=amount/6.02214199E+26;

   % Assign dummy res and rnfwhm
   res=zeros(nchan,1);
   rnfwhm=zeros(nchan,1);


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Re-order the ctrans layers to run top to bottom
   if (ip ==1)
      indr=nlay:-1:1;
      nanglaygas=nang*nlay*ngas;
      nanglay=nang*nlay;
      indis=1:nanglaygas;
      indwant=indis;
      ilay=1:nlay;
      ileft=nang*(ilay-1);
      iright=nang*(nlay - ilay);
      for igas=1:ngas
         ioffset=(igas - 1)*nanglay;
         for ia=1:nang
            indwant(ioffset + ileft + ia)=indis(ioffset + iright + ia);
         end
      end
   end
   roctrans=ctrans(:,indwant)'; %'

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Write out the merged profile + transmittance data file
   outname=[outfile int2str(ipnum) '.dat'];
   titlecom=[comment ' p#' int2str(ipnum)];
if (lfcow2fow == 1 )
[iok]=wrt_convdat_fcow2fow(outname, nang, nlay, ngas, nchan, gasids, secang,...
titlecom, temp, amount, fchan, ichan, res, rnfwhm, roctrans);
else
   [iok]=wrt_convdat(outname, nang, nlay, ngas, nchan, gasids, secang, ...
      titlecom, temp, amount, fchan, ichan, res, rnfwhm, roctrans);
end
   if (iok == 0)
      error(['Error detected writing ' outname]);
   end

end
ok=1;

%%% end of function %%%
