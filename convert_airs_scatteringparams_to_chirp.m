%% See ~/asl/sergiodisk/WorkDirDec2025/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/srcF77_jac/fnmie_iceGHMbaum_waterdrop_desertdust_chirp.f

boo = load('/home/sergio/asl/rta/sarta_database/Data_mie_apr08/chirp_ichan_vchan.mat');
fchirp = boo.h2.vchan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ICE

infile =  '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/ice_baumGHM_abs.dat';
outfile = '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/ice_baumGHM_abs_chirp.dat';
interp_sarta_scatter_coeffs_old_instr_to_new_instr(infile,fchirp,outfile);

infile =  '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/ice_baumGHM_ext.dat';
outfile = '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/ice_baumGHM_ext_chirp.dat';
interp_sarta_scatter_coeffs_old_instr_to_new_instr(infile,fchirp,outfile);

infile =  '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/ice_baumGHM_asy.dat';
outfile = '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/ice_baumGHM_asy_chirp.dat';
interp_sarta_scatter_coeffs_old_instr_to_new_instr(infile,fchirp,outfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% water

infile =  '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/waterdrop_abs.dat';
outfile = '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/waterdrop_abs_chirp.dat';
interp_sarta_scatter_coeffs_old_instr_to_new_instr(infile,fchirp,outfile);

infile =  '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/waterdrop_ext.dat';
outfile = '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/waterdrop_ext_chirp.dat';
interp_sarta_scatter_coeffs_old_instr_to_new_instr(infile,fchirp,outfile);

infile =  '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/waterdrop_asy.dat';
outfile = '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/waterdrop_asy_chirp.dat';
interp_sarta_scatter_coeffs_old_instr_to_new_instr(infile,fchirp,outfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% dust

infile =  '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/desertdust_abs.dat';
outfile = '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/desertdust_abs_chirp.dat';
interp_sarta_scatter_coeffs_old_instr_to_new_instr(infile,fchirp,outfile);

infile =  '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/desertdust_ext.dat';
outfile = '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/desertdust_ext_chirp.dat';
interp_sarta_scatter_coeffs_old_instr_to_new_instr(infile,fchirp,outfile);

infile =  '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/desertdust_asy.dat';
outfile = '/home/sergio/asl/rta/sarta_database/Data_mie_apr08/desertdust_asy_chirp.dat';
interp_sarta_scatter_coeffs_old_instr_to_new_instr(infile,fchirp,outfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str = ['ls -lt /home/sergio/asl/rta/sarta_database/Data_mie_apr08/*chirp.dat'];
eval(str)

