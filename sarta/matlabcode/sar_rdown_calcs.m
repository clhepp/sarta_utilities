% 

cd /home/chepplew/projects/sarta/cris_hr

addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil

 [hdr har pdr par] = rtpread('regr49_1100_400ppm_2235g4.op.rtp');

[nx ny] = size(pdr.plevs)
pdr.upwell = single(2.*ones(1,ny));
pdr.pobs   = single(1013.25*ones(1,ny));
pdr.zobs   = single(ones(1,ny));

 fortp='regr49_1100_400_2235g4_rdown.rtp';
 rtpwrite(fortp,hdr,har,pdr,par);

%{
% now run sarta: 
n72> SARTA_EXE=/home/chepplew/gitLib/sarta/bin/sarta_kc_cris_hrg4_400p_full_th2
n72> FIN=/home/chepplew/projects/sarta/cris_hr/regr49_1100_400_2235g4_rdown.rtp
n72> FOUT=/asl/s1/chepplew/projects/sarta/cris_hr/sar_r49_rdown.rtp
n72> rm sarta.out
n72> $SARTA_EXE fin=$FIN fout=$FOUT 1>sarta.out
%}

% see what the result look like, compare with a Scott file
fnsar = '/asl/s1/chepplew/projects/sarta/cris_hr/sar_r49_rdown.rtp';
[hds has pds pas] = rtpread(fnsar);
 bsc = real(rad2bt(hds.vchan,pds.rcalc));

