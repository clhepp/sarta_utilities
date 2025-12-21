function [fIASI,dme,aeg] = read_scatter_coeffs(fileIASI)

%% see eg SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/srcF77_jac/rdcldt.f

%% C
%% C         Read mie particle sizes
%%           READ(IOUN) (MIEPS(IL,K),IL=1,MIENPS(K))
%% c          print *,'READ MIESIZES'
%% C
%% C         Read mie abs data for required channels
%%           DO I=1,MXCHAN
%%              IF (INDCHN(I) .NE. 0) THEN
%% 		READ(IOUN) (MIEABS(INDCHN(I),IL,K),IL=1,MIENPS(K))
%%              ELSE
%%                 READ(IOUN) (XJUNK(IL),IL=1,MIENPS(K))
%%              ENDIF
%%           ENDDO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(fileIASI,'r','ieee-be');

%% read nchan and ndme
irecsize = fread(fid,1,'integer*4');
fprintf(1,'N irecsize = %5i \n',irecsize);
nchan    = fread(fid,1,'integer*4');
ndme     = fread(fid,1,'integer*4');
irecsize = fread(fid,1,'integer*4');
fprintf(1,'N irecsize = %5i \n',irecsize);

iFound = 0;
if nchan ~= 8461
  iFound = iFound - 1;
  %disp('need nchan = 8461 for IASI')
else
  disp(' file has 8461 channels ==> IASI')
  iFound = iFound + 1;
  %% https://space.oscar.wmo.int/instruments/view/iasi
  iX = 1 : 8461;
  dv  = 0.25;
  fIASI = 645 + (iX-1)*dv;
end

if nchan ~= 2834
  iFound = iFound - 1;
  %disp('need nchan = 2834 for AIRS')
else
  disp(' file has 2834 channels ==> AIRS')
  iFound = iFound + 1;
  boo = load('/home/sergio/asl/rta/sarta_database/Data_mie_apr08/chanid_freq_m140.txt');
  fIASI = boo(:,2);
end

if nchan ~= 1702
  iFound = iFound - 1;
  %disp('need nchan = 1702 for CHIRP')
else
  disp(' file has 1702 channels ==> chirp')
  iFound = iFound + 1;
  boo = load('/home/sergio/asl/rta/sarta_database/Data_mie_apr08/chirp_ichan_vchan.mat');
  fIASI = boo.h2.vchan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% read dme
irecsize=fread(fid,1,'integer*4');
fprintf(1,'D irecsize = %5i \n',irecsize)
dme = fread(fid,ndme,'real*4');
irecsize=fread(fid,1,'integer*4');
fprintf(1,'D irecsize = %5i \n',irecsize)
fprintf(1,'ndme = %3i \n',ndme)

%% read the mie coeff (could be aeg == absorption/extinction/asymmetry)
for ii = 1 : nchan
  irecsize=fread(fid,1,'integer*4');
  if ii == 1
    fprintf(1,'X irecsize = %5i \n',irecsize)
  end
  junk = fread(fid,ndme,'real*4');
  irecsize=fread(fid,1,'integer*4');
  if ii == 1  
    fprintf(1,'X irecsize = %5i \n',irecsize)
  end
  aeg(ii,:) = junk;
end

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pcolor(fIASI,dme,aeg'); shading interp; colorbar; pause(0.1);
