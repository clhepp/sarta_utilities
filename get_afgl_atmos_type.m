function [iiATM] = get_afgl_atmos_type(profin)
%
% Get index (1 to 6) for atmosphere type of profin
% with afgl types: TRP,MLS,MLW,SAS,SAW,STD (1 to 6)
%

addpath /asl/matlib/time

[yy,mon,day,hr,~,~] = datevec(tai2dtime(profin.rtime));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iaTRP = find(abs(profin.rlat) < 30);

iaMLS_N = find(profin.rlat < +60  & profin.rlat >= +30 & mon >= 6 & mon <= 9);
iaMLS_S = find(profin.rlat > -60  & profin.rlat <= -30 & mon >= 1 & mon <= 4);
iaMLS = union(iaMLS_N,iaMLS_S);

iaMLW_N = find(profin.rlat < +60  & profin.rlat >= +30 & mon >= 1 & mon <= 4);
iaMLW_S = find(profin.rlat > -60  & profin.rlat <= -30 & mon >= 6 & mon <= 9);
iaMLW = union(iaMLW_N,iaMLW_S);

iaSTD = find(abs(profin.rlat) < 60 & abs(profin.rlat) >= 30);
iaSTD = setdiff(setdiff(iaSTD,iaMLS),iaMLW);

iaSAS_N = find(profin.rlat >= +60 & mon >= 4 & mon <= 9);
iaSAS_S = find(profin.rlat <= -60 & ((mon >= 10 & mon <= 12) | (mon >= 1 & mon <= 3)));
iaSAS = union(iaSAS_N,iaSAS_S);

iaSAW_N = find(profin.rlat >= +60 & ((mon >= 10 & mon <= 12) | (mon >= 1 & mon <= 3)));
iaSAW_S = find(profin.rlat <= -60 & mon >= 4 & mon <= 9);
iaSAW = union(iaSAW_N,iaSAW_S);

% convert to group number:
iiATM = zeros(size(profin.rlat));
iiATM(iaTRP) = 1;
iiATM(iaMLS) = 2;
iiATM(iaMLW) = 3;
iiATM(iaSAS) = 4;
iiATM(iaSAW) = 5;
iiATM(iaSTD) = 6;

