% compare_thermal_coefs.m

addpath /home/chepplew/projects/sarta/matlabcode       % rdcoef.m

FNCOF1='/asl/data/sarta_database/Data_CrIS_apr09/Coef/therm.dat';
[ich1 fch1 coef1 info1] = rdcoef(11,0,FNCOF1);

FNCOF2='/asl/s1/chepplew/data/sarta_database/Data_CrIS_hrg4_400ppm/Coef/therm_05082016.dat';
[ich2 fch2 coef2 info2] = rdcoef(11,0,FNCOF2);


FNCOF3='/asl/data/sarta_database/Data_AIRS_apr08/Coef/therm_m140.dat';
[ich3 fch3 coef3 info3] = rdcoef(11,0,FNCOF3); 


figure(4);clf;hold on; for i=1:6 plot(fch1,coef1(:,1,i),'-'); end, grid on;axis([640 1100 -3 3]);
figure(5);clf;hold on; for i=1:6 plot(fch2,coef2(:,1,i),'-');end;grid on;axis([640 1100 -Inf Inf]);

