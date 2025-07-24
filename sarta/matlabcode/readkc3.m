function [data, wnums] = readkc3(kfile, dfile)

% function [data, wnums] = readkc3(kfile, dfile);
%
% readkc3 is a simple reader & unchunker for kcarta output files.
% Identical to readkcstd except array dimensions nrow & ncol (each a
% 4 byte integer) have been added to the start of dfile.
%
% INPUTS:
%    kfile  - name of kcarta output file
%    dfile  - name of simple binary output file to create (optional)
%
% OUTPUTS:
%    data   - (w x n) array of data from kcarta
%    wnums  - (w x 1) vector of data wavenumbers
%
% If the input parameter dfile is specified, then the data array
% is written to file dfile, and the return values [data, wnums]
% are left unassigned.
%
% The w by n data array has one row for each output wavenumber and
% one column for each "output data block" row.  In practice, this
% means the data array columns are a concatenation of whatever was 
% asked for in the kcarta driver file.
% 
% LIMITATIONS
% 
% Large data sets (e.g., mixed path sets of more than one or two
% "chunks") take up too much space to be returned in memory; for 
% such cases an output file should always be specified.

% H. Motteler, 8/24/98
% Updated for version 1.03 by Scott Hannon, 30 July 1999
% Updated by Sergio (when?) to also handle long header data files
% 
% Scott Hannon, 21 July 2000, for v1.06, add nrow & ncol to start of dfile
% Version readkc3.m created 8 April 2002, by Scott Hannon.  Based on
%    readkc2.m, but with freq points added to output file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fin,msg] = fopen(kfile, 'r');  % Open in native format
if fin == -1
  error(['error opening input file\n', msg]);
  end
fid=fin; % Sergio's modification

%%%%%%%%%%%%%%%%%%%%%%
% READ HEADER RECORDS
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%
% MAIN HEADER

% version number
flen    = fread(fin, 1, 'integer*4');
version = fread(fin, 80, 'char');
version = setstr(version');
flen    = fread(fin, 1, 'integer*4');

% number of layers
flen   = fread(fin, 1, 'integer*4');
nlayer = fread(fin, 1, 'integer*4');
flen   = fread(fin, 1, 'integer*4');
iNumLayers=nlayer;  % Sergio's modification

% number of params 
flen    = fread(fin, 1, 'integer*4');
nparams = fread(fin, 1, 'integer*4');
flen    = fread(fin, 1, 'integer*4');
flen    = fread(fin, 1, 'integer*4');
rparams = fread(fin, nparams, 'real*4');
flen    = fread(fin, 1, 'integer*4');

% comment
flen    = fread(fin, 1, 'integer*4');
comment = fread(fin, 80, 'char');
comment = setstr(comment');
flen    = fread(fin, 1, 'integer*4');

% start, stop frequency
flen = fread(fin, 1, 'integer*4');
fmin = fread(fin, 1, 'real*4');
fmax = fread(fin, 1, 'real*4');
flen = fread(fin, 1, 'integer*4');

% low and high chunk index
flen      = fread(fin, 1, 'integer*4');
lowchunk  = fread(fin, 1, 'integer*4');
highchunk = fread(fin, 1, 'integer*4');
flen      = fread(fin, 1, 'integer*4');

% read iLongOrShort
flen  = fread(fin, 1, 'integer*4');
htype = fread(fin, 1, 'integer*4');
flen  = fread(fin, 1, 'integer*4');

%if htype ~= -1
%  error('readkc: can only do short headers');
%  end

%%%%
% This block added by Scott for version 1.03+
version_number=str2num(version(2:5));
if ((version_number >= 1.03) & (version_number <= 1.05)) 
  % read in M1000mb,M100mb,MSubLayer,M50mb,M10mb,MThickLayer
  flen = fread(fin,1,'integer*4');
  iaS  = fread(fin,6,'integer*4');
  flen = fread(fin,1,'integer*4');

  % read in pressure levels
  flen = fread(fin,1,'integer*4');
  raP  = fread(fin,nlayer+1,'real*4'); % Note: nlayer+1
  flen = fread(fin,1,'integer*4');

% Larry McMillin found bug in pProf writing out one too many elements
elseif (version_number >= 1.06)
  % read in M1000mb,M100mb,MSubLayer,M50mb,M10mb,MThickLayer
  flen = fread(fin,1,'integer*4');
  iaS  = fread(fin,6,'integer*4');
  flen = fread(fin,1,'integer*4');

  % read in pressure levels
  flen = fread(fin,1,'integer*4');
  raP  = fread(fin,nlayer,'real*4');   % Note: nlayer
  flen = fread(fin,1,'integer*4');
  end


%%%%%%%%%%%%%%%%%
% GAS PATH HEADER

if htype < 0 %%% short version %%%
  % number of paths to be output
  flen    = fread(fin, 1, 'integer*4');
  ngasout = fread(fin, 1, 'integer*4');
  flen    = fread(fin, 1, 'integer*4');

  if (ngasout > 0)
    flen     = fread(fin, 1, 'integer*4');
    gaspaths = fread(fin, ngasout, 'integer*4');
    flen     = fread(fin, 1, 'integer*4');
    end
  end

if htype > 0 %%% long version ..... copied direct from readgaspaths.m %%%
  %Note: assumes iNumPaths=iNumGases*iNumLayers 
  flen       = fread(fid,1,'integer*4'); 
  iNumPaths  = fread(fid,1,'integer*4'); 
  flen       = fread(fid,1,'integer*4'); 
  ngasout    = iNumPaths;  % Sergio's modification
   
  iNumGases = iNumPaths/iNumLayers;   % This is number of gases 
  iaGasID   = zeros(1,iNumGases); 
  raaAmt    = zeros(iNumLayers,iNumGases); 
  raaTemp   = zeros(iNumLayers,iNumGases); 
 
  % Now read the PathNum GasID Temperature Amount for each path in the profile 
  for ii=1:iNumGases 
    for jj=1:iNumLayers 
      flen   = fread(fid,1,'integer*4'); 
      iPath  = fread(fid,1,'integer*4'); 
      iGasID = fread(fid,1,'integer*4'); 
      rTemp  = fread(fid,1,'real*4'); 
      rAmt   = fread(fid,1,'real*4'); 
      flen   = fread(fid,1,'integer*4'); 
      raaAmt(jj,ii)  = rAmt; 
      raaTemp(jj,ii) = rTemp; 
      end  
    iaGasID(ii)=iGasID; 
    end 
 
  % Read in Number of Paths to be output each time 
  flen         = fread(fid,1,'integer*4'); 
  iNumPathsOut = fread(fid,1,'integer*4'); 
  flen         = fread(fid,1,'integer*4'); 
  if ((iNumPathsOut < 0) | (iNumPathsOut > iNumPaths)) 
    error('Mistake in number of paths to be output!!!') 
    end 
  ia=0; 
  clear ia; 
  if (iNumPathsOut > 0) 
    % Read in paths 
    iaOutPaths = zeros(1,iNumPathsOut); 
    flen       = fread(fid,1,'integer*4'); 
    ia         = fread(fid,iNumPathsOut,'integer*4'); 
    iaOutPaths(1,:)=ia'; 
    clear ia 
    flen       = fread(fid,1,'integer*4'); 
    gaspaths   = iaOutPaths; % Sergio's modification
    end 
  end %%% end of long version htype > 0 %%%


%%%%%%%%%%%%%%%%%%%
% MIXED PATH HEADER
if htype < 0 %%% short version %%%
  flen      = fread(fin, 1, 'integer*4');
  nmixpaths = fread(fin, 1, 'integer*4');
  flen      = fread(fin, 1, 'integer*4');

  if (nmixpaths > 0)

    flen    = fread(fin, 1, 'integer*4');
    nmixout = fread(fin, 1, 'integer*4');
    flen    = fread(fin, 1, 'integer*4');

    if (nmixout > 0)
      flen     = fread(fin, 1, 'integer*4');
      mixpaths = fread(fin, nmixout, 'integer*4');
      flen     = fread(fin, 1, 'integer*4');
      end
    end
  end

if htype > 0 %%% long version ... copied direct from readmixedpaths.m %%%
  % Read in number of mixed paths in *MIXFIL 
  flen      = fread(fid,1,'integer*4'); 
  iNpmix    = fread(fid,1,'integer*4'); 
  flen      = fread(fid,1,'integer*4'); 
  nmixpaths = iNpmix; % Sergio's modification

  % Read in #of lines that have info in *MIXFIL 
  flen           = fread(fid,1,'integer*4'); 
  iMixFileLines  = fread(fid,1,'integer*4'); 
  flen           = fread(fid,1,'integer*4'); 
 
  if (iNpmix > 0) 
    caaMixPathInfo = zeros(iMixFileLines,130); 
    % Read in the mixing table 
    for ii=1:iMixFileLines 
      flen           = fread(fid,1,'integer*4'); 
      caMixPathInfo  = fread(fid,130,'char'); 
      flen           = fread(fid,1,'integer*4'); 
      end 
 
    % Read in mixed path temperatures 
    raMixTemp = zeros(1,iNpmix); 
    flen      = fread(fid,1,'integer*4'); 
    ra        = fread(fid,iNpmix,'real*4'); 
    flen      = fread(fid,1,'integer*4'); 
    raMixTemp(1,:)=ra(1:length(ra))'; 
    clear ra  
 
    % Read in Number of mixed paths to be output each time 
    flen            = fread(fid,1,'integer*4'); 
    iNumMixPathsOut = fread(fid,1,'integer*4'); 
    flen            = fread(fid,1,'integer*4'); 
    nmixout         = iNumMixPathsOut; % Sergio's modification

    if ((iNumMixPathsOut < 0) | (iNumMixPathsOut > iNpmix)) 
      error('Mistake in number of mixed paths to be output!!!') 
      end 
    iaOutMixPaths=0; 
    clear iaOutMixPaths; 
    if (iNumMixPathsOut > 0) 
      % Read in paths 
      iaOutMixPaths = zeros(1,iNumMixPathsOut); 
      flen          = fread(fid,1,'integer*4'); 
      ia            = fread(fid,iNumMixPathsOut,'integer*4'); 
      flen          = fread(fid,1,'integer*4'); 
      %fprintf(1,'%3i\n',iaOutMixPaths); 
      iaOutMixPaths(1,:)=ia'; 
      clear ia   
      mixpaths      = iaOutMixPaths; % Sergio's modification
      end 
    end 
  end %%% end of long version LS > 0 %%%


%%%%%%%%%%%%%%%%%%%
% ATMOSPHERE HEADER

% number of atmospheres in *RADFIL
flen   = fread(fin, 1, 'integer*4');
natmos = fread(fin, 1, 'integer*4');
flen   = fread(fin, 1, 'integer*4');

% max number of emissivity points per atm
flen  = fread(fin, 1, 'integer*4');
nemis = fread(fin, 1, 'integer*4');
flen  = fread(fin, 1, 'integer*4');

if (natmos > 0)

  for i = 1:natmos

    % read in atmosphere number, # of mixed paths in atmosphere
    flen = fread(fin, 1, 'integer*4');
    iAtm = fread(fin, 1, 'integer*4');
    iNumPathsAtmos = fread(fin, 1, 'integer*4');
    flen = fread(fin, 1, 'integer*4');

    if (iNumPathsAtmos > 0)

      % Read in paths making up atmosphere
      flen = fread(fin, 1, 'integer*4');
      ia   = fread(fin, iNumPathsAtmos, 'integer*4');
      flen = fread(fin, 1, 'integer*4');

      flen   = fread(fin, 1, 'integer*4');
      rTbdy  = fread(fin, 1, 'real*4');
      rTinit = fread(fin, 1, 'real*4');
      rSat   = fread(fin, 1, 'real*4');
      rHgt   = fread(fin, 1, 'real*4');
      flen   = fread(fin, 1, 'integer*4');

      flen           = fread(fin, 1, 'integer*4');
      ikSolar        = fread(fin, 1, 'integer*4');
      rkSolarAngle   = fread(fin, 1, 'real*4');
      rkSolarRefl    = fread(fin, 1, 'real*4');
      ikThermal      = fread(fin, 1, 'integer*4');
      rkThermalAngle = fread(fin, 1, 'real*4');
      ikThermalJacob = fread(fin, 1, 'integer*4');
      flen           = fread(fin, 1, 'integer*4');

      flen = fread(fin, 1, 'integer*4');
      iEms = fread(fin, 1, 'integer*4');
      flen = fread(fin, 1, 'integer*4');

      % Now read the freq, ems data points
      for jj=1:iEms
        flen = fread(fin, 1, 'integer*4');
        rff  = fread(fin, 1, 'real*4');
        ree  = fread(fin, 1, 'real*4');
        flen = fread(fin, 1, 'integer*4');
        end 

      % Now read in layers to be output
      flen          = fread(fin, 1, 'integer*4');
      iNumLayersOut = fread(fin, 1, 'integer*4');
      flen          = fread(fin, 1, 'integer*4');

      if (iNumLayersOut > 0)
        % read in mixed paths to be output
        flen = fread(fin, 1, 'integer*4');
        ia   = fread(fin, iNumLayersOut, 'integer*4');
        flen = fread(fin, 1, 'integer*4');
        end

      if (iNumLayersOut > 0)
        % read in pressures at which radiances to be output
        flen = fread(fin, 1, 'integer*4');
        ra   = fread(fin, iNumLayersOut, 'real*4');
        flen = fread(fin, 1, 'integer*4');
        end
      end
    end
  end

%%%%%%%%%%%%%%%%%%%%%%
% DATA RECORD SUMMARY
%%%%%%%%%%%%%%%%%%%%%%

flen   = fread(fin, 1, 'integer*4');
nchunk = fread(fin, 1, 'integer*4'); % total number of chunks
nODBs  = fread(fin, 1, 'integer*4'); % number of (ODBs) "output data blocks"
flen   = fread(fin, 1, 'integer*4');

flen1    = fread(fin, 1, 'integer*4');
nODBrows = fread(fin, nODBs, 'integer*4'); % rows in each ODB
flen2    = fread(fin, 1, 'integer*4');

% sanity checks
if flen1 ~= flen2
  error('Fortran records out of phase!');
end

if nchunk ~= 1+highchunk-lowchunk
  error('readkc: chunk counts do not match!');
end

fprintf(2, 'readkc: %d chunks, %d ODBs, %d total rows\n', ...
        nchunk, nODBs, sum(nODBrows));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize output file or array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noutrow = nchunk * 10000;       % number of output rows
noutcol = sum(nODBrows);        % number of output columns
ODBrowdat = zeros(10000, 1);    % ODB row data buffer

if nargin == 1
 
  % initialize output arrays
  data = zeros(noutrow, noutcol);
  wnums = zeros(noutrow, 1);

else

  % pre-extend the output file
%  totbytes = noutrow * noutcol * 4;  % never used so not update
  fprintf(2, 'readkc: creating %d x %d element output array\n', ...
          noutrow, noutcol);

  [fout,msg] = fopen(dfile, 'w');
  if fout == -1
    error(['error creating output file', msg]);
  end

  %%%
  % Scott's mod to add nrow & ncol to start of dfile
  % Output mini-header = noutrow and noutcol
  fwrite(fout, noutrow, 'integer*4');
  fwrite(fout, noutcol, 'integer*4');
  %%%

  z = zeros(noutrow,1);

  %%%
  % Scott's mod for freq points
  if fwrite(fout, z, 'real*8') ~= noutrow
     error('fwrite in extend failed for freq points');
  end
  %%%

  for i = 1:noutcol
    if fwrite(fout, z, 'real*4') ~= noutrow
      error(['fwrite in extend failed, i=',num2str(i)]);
    end
  end
  fclose(fout);
  clear z

  % open the extended file in "update" mode
  [fout,msg] = fopen(dfile, 'r+');
  if fout == -1
    error(['error opening output file', msg]);
  end
  disp('finished pre-extending output file')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read "output data blocks"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% Variables needed to shift file position on account of nrow,ncol,
% and freq data added to start of file by Scott
offset1=4 + 4;                     % nrow and ncol are both integer*4
offset2=offset1 + nchunk*10000*8;  % freq points are real*8
%%%

for chunk = 1:nchunk
disp(['working on chunk ' int2str(chunk)])
  cumODBrow = 1;

  for odb = 1:nODBs

    % read ODB header
    flen    = fread(fin, 1, 'integer*4');
    type    = fread(fin, 1, 'integer*4'); % ODB type
    subtype = fread(fin, 1, 'integer*4'); % ODB subtype
    nrow    = fread(fin, 1, 'integer*4'); % number of rows, this ODB
    flen    = fread(fin, 1, 'integer*4');

    flen    = fread(fin, 1, 'integer*4');
    ncols   = fread(fin, 1, 'integer*4'); % spectral pts in a chunk
    frlow   = fread(fin, 1, 'real*4');    % low freq of chunk
    frhigh  = fread(fin, 1, 'real*4');    % high freq of chunk
    frinc   = fread(fin, 1, 'real*4');    % frequency increment
    flen    = fread(fin, 1, 'integer*4');
  
    if chunk == 1
      fprintf(2, 'readkc: ODB type = %d, subtype = %d, rows = %d\n', ...
                type, subtype, nrow);
    end

    % Freq info; Scott
    flo=round(frlow*100)/100;  % round to nearest 0.01 cm-1
    df=round( (frhigh-frlow)*100 )/1000000;  % 100*10000

    % sanity check
    if nrow ~= nODBrows(odb)
      fprintf(2, 'readkc: WARNING -- bad ODB row count,  ');
      fprintf(2, 'header = %-4d ODB = %-4d\n', nODBrows(odb), nrow);
    end

    if (nargin == 2 & odb == 1)
       % write frequency points for current chunk
       outpos = offset1 + (chunk-1)*10000*8;
%       wnums = frlow + (0:9999)*frinc; % not accurate enough
       wnums = flo + (0:9999)*df;
       if fseek(fout, outpos, -1) ~= 0
          error(['fseek failed, odb=',num2str(odb),' chunk=',num2str(chunk)]);
       end
       if fwrite(fout, wnums, 'real*8') ~= 10000
          error(['fwrite failed, odb=',num2str(odb),' chunk=',num2str(chunk)]);
       end
    end


    % read a row of ODB data
    for ODBrow = 1:nODBrows(odb)
      flen1 = fread(fin, 1, 'integer*4');
      [ODBrowdat, count] = fread(fin, [10000,1], 'real*4');
      if count ~= 10000
        error(['fread failed, odb=',num2str(odb),' chunk=',num2str(chunk)]);
      end
      flen2 = fread(fin, 1, 'integer*4');

      if nargin == 1

        % write ODBrowdat to the data array
        outrows = (chunk-1) * 10000 + 1 : chunk * 10000;
        outcol = cumODBrow;
        data(outrows, outcol) = ODBrowdat;
%        wnums(outrows) = frlow + (0:9999)*frinc;
        wnums(outrows) = flo + (0:9999)*df;  % Scott
      else

        % write ODBrowdat to the output file
        % modified by Scott for nrow, ncol, wnums

        %outpos = ((cumODBrow-1)*nchunk*10000 + (chunk-1)*10000)*4;
        %outpos = 8 + ((cumODBrow-1)*nchunk*10000 + (chunk-1)*10000)*4;
        outpos = offset2 + ((cumODBrow-1)*nchunk*10000 + (chunk-1)*10000)*4;
        if fseek(fout, outpos, -1) ~= 0
          error(['fseek failed, odb=',num2str(odb),' chunk=',num2str(chunk)]);
        end
        if fwrite(fout, ODBrowdat, 'real*4') ~= 10000
          error(['fwrite failed, odb=',num2str(odb),' chunk=',num2str(chunk)]);
        end     

      end

      cumODBrow = cumODBrow + 1;

    end % loop on current ODB rows

  end % loop on ODBs

end % loop on chunks

fclose(fin);

if nargin == 2 
  fclose (fout);
end

%%% end of function %%%
