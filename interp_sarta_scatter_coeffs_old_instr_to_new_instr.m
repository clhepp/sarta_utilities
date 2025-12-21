function interp_sarta_scatter_coeffs_old_instr_to_new_instr(fileIASI,freqNEW,fileNEW)

%% takes in fileIAS1 : string  (since it will have higest wavenumber resolution)
%% reads in contents
%% interpolats to wavneumber freqNEW
%% saves to fileNEW
%{
boo = load('/home/sergio/asl/rta/sarta_database/Data_mie_apr08/chirp_ichan_vchan.mat');
fchirp = boo.h2.vchan;
infile =  '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/ice_baumGHM_abs.dat';
outfile = '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/ice_baumGHM_abs_chirp.dat';
interp_sarta_scatter_coeffs_old_instr_to_new_instr(infile,fchirp,outfile);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(fileNEW)
  str = [fileNEW ' already exists'];
  error(str);
end

if ~exist(fileIASI)
  str = [fileIASI ' does not exist'];
  error(str);
end

freqNEW0 = freqNEW;
freqNEW = unique(freqNEW);
good = find(isfinite(freqNEW));
freqNEW = freqNEW(good);

if length(freqNEW0) ~= length(freqNEW)
  [length(freqNEW0) length(freqNEW)]
  disp('length(freqNEW0) ~ length(freqNEW)')
end

[fIASI,dme,aeg] = read_sarta_scatter_coeffs(fileIASI);
for ii = 1 : length(dme)
  aegNew(:,ii) = interp1(fIASI,aeg(:,ii),freqNEW,[],'extrap');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);
  pcolor(freqNEW,dme,aegNew'); shading interp; colorbar; pause(0.1);

whos fIASI dme aeg freqNEW dme aegNew

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(fileNEW,'w','ieee-be');

nchan = length(freqNEW);
ndme = length(dme);

fwrite(fid,4*2,'integer*4');
fwrite(fid,nchan,'integer*4');
fwrite(fid,ndme,'integer*4');
fwrite(fid,4*2,'integer*4');

fwrite(fid,4*ndme,'integer*4');
fwrite(fid,dme,'real*4');
fwrite(fid,4*ndme,'integer*4');

for ii = 1 : nchan
  junk = aegNew(ii,:);
  fwrite(fid,4*ndme,'integer*4');
  fwrite(fid,junk,'real*4');
  fwrite(fid,4*ndme,'integer*4');  
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iCheck = -1;
if iCheck > 0
  disp('checking everything is OK')
  [fx,dx,aegx] = read_sarta_scatter_coeffs(fileNEW);
  whos freqNEW dme aegNew
  whos fx dx aegx
  [sum(freqNEW-fx) sum(dme-dx) sum(aegNew(:)-aegx(:))]
  figure(3)
    pcolor(freqNEW,dme,aegNew'-aegx'); shading interp; colorbar; pause(0.1);
end
