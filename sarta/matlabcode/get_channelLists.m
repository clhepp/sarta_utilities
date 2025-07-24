% get_channelLists.m

% CLH May 2016
% Extract CrIS hi_res (080808) channel lists for use with SARTA fit_ftc development
% with 4 guard channels

% cd /home/chepplew/gitLib/chanLists

% First the full bands
B1 = [1:721; 647.500:0.625:1097.5];
B2 = [722:1594; 1207.500:0.625:1752.500];
B3 = [1595:2235; 2152.500:0.625:2552.5];  

opFn = fopen('list_crisB1_hrg4','w');
  fprintf(opFn, '%4i\t %8.3f\n',B1);
fclose(opFn)
opFn = fopen('list_crisB2_hrg4','w') ;
   fprintf(opFn, '%4i\t %8.3f\n',B2);
fclose(opFn);
opFn = fopen('list_crisB3_hrg4','w');
   fprintf(opFn, '%4i\t %8.3f\n',B3);
fclose(opFn);

% Now the Sets - derived from Scott's sets of coefficients which he defined 
% for 2 guard channels per edge.

addpath /home/chepplew/projects/sarta/matlabcode/
fname1='/asl/s1/chepplew/data/sarta_database/Data_CrIS_HR_mar16/Coef/set1.dat';
   [ichan1, fchan1, coef1, info1] = rdcoef(1,0, fname);
fname2='/asl/s1/chepplew/data/sarta_database/Data_CrIS_HR_mar16/Coef/set2_hrg2.dat';
   [ichan2, fchan2, coef2, info2] = rdcoef(1,0, fname2);   
fname3='/asl/s1/chepplew/data/sarta_database/Data_CrIS_HR_mar16/Coef/set3_hrg2.dat';
   [ichan3, fchan3, coef3, info3] = rdcoef(3,0, fname3);  
fname4='/asl/s1/chepplew/data/sarta_database/Data_CrIS_HR_mar16/Coef/set4_hrg2.dat';
   [ichan4, fchan4, coef4, info4] = rdcoef(4,0, fname4); 
fname5='/asl/s1/chepplew/data/sarta_database/Data_CrIS_HR_mar16/Coef/set5_hrg2.dat';
   [ichan5, fchan5, coef5, info5] = rdcoef(5,0, fname5);
fname6='/asl/s1/chepplew/data/sarta_database/Data_CrIS_HR_mar16/Coef/set6_hrg2.dat';
   [ichan6, fchan6, coef6, info6] = rdcoef(6,0, fname6);        
fname7='/asl/s1/chepplew/data/sarta_database/Data_CrIS_HR_mar16/Coef/set7_hrg2.dat';
   [ichan7, fchan7, coef7, info7] = rdcoef(7,0, fname7);                               


junk = ichan1+2;                                                                   
ichan1g4=[1 2 junk']';                                                             
fchan1g4 = [647.500 648.125 fchan1']';

junk = ichan2+2;                                                                   
ichan2g4 = [junk' 720 721]';                               
fchan2g4 = [fchan2' 1096.875 1097.500]';

junk = ichan3+6;
ichan3g4 = [722 723 junk']';
junk = [ichan3g4' 1593 1594]';
ichan3g4 = junk;
fchan3g4 = [1207.500 1208.125 fchan3' 1751.875 1752.500]';

junk = ichan4+10;             
ichan4g4 = [1595 1596 junk']';
fchan4g4 = [2152.500 2153.125 fchan4']';

ijunk = ichan5+10;
ichan5g4 = [ijunk' 2234 2235]';
fchan5g4 = [fchan5' 2551.875 2552.500]';

ichan6g4 = ichan6+10;

ichan7g4 = ichan7+10;

opFn = fopen('list_crisSet1_hrg4','w')
  for i=1:493; fprintf(opFn, '%4i\t %8.3f\n',ichan1g4(i),fchan1g4(i)); end
fclose(opFn)
opFn = fopen('list_crisSet2_hrg4','w')                                         
  for i=1:228; fprintf(opFn, '%4i\t %8.3f\n',ichan2g4(i),fchan2g4(i)); end
fclose(opFn)                                                            
opFn = fopen('list_crisSet3_hrg4','w');                                   
  for i=1:873; fprintf(opFn, '%4i\t %8.3f\n',ichan3g4(i),fchan3g4(i)); end
fclose(opFn)
opFn = fopen('list_crisSet4_hrg4','w');                                                  
  for i=1:126; fprintf(opFn, '%4i\t %8.3f\n',ichan4g4(i),fchan4g4(i)); end               
fclose(opFn);
opFn = fopen('list_crisSet5_hrg4','w');                                                  
   for i=1:75; fprintf(opFn, '%4i\t %8.3f\n',ichan5g4(i),fchan5g4(i)); end                
fclose(opFn)
opFn = fopen('list_crisSet6_hrg4','w');                                  
   for i=1:415; fprintf(opFn, '%4i\t %8.3f\n',ichan6g4(i),fchan6(i)); end  
fclose(opFn);
opFn = fopen('list_crisSet7_hrg4','w');                                   
   for i=1:25; fprintf(opFn, '%4i\t %8.3f\n',ichan7g4(i),fchan7(i)); end   
fclose(opFn);                                                        

%{
figure(2);clf;plot(ichan1g4,fchan1g4,'b+');grid on;hold on;
plot(ichan2g4,fchan2g4,'r+');                              

figure(4);clf;plot(ichan3g4, fchan3g4,'o');grid on;

figure(5);clf;plot(ichan4g4,fchan4g4,'bo');grid on;axis([1590 2240 2150 2560]);hold on;
 plot(ichan5g4,fchan5g4,'go');                                                          
 plot(ichan6g4,fchan6,'ro');                                                            
 plot(ichan7g4,fchan7,'mo');                                                            

 L1=textread('list_crisSet1_hrg4');    
 L2=textread('list_crisSet2_hrg4');
 L3=textread('list_crisSet3_hrg4');
 L4=textread('list_crisSet4_hrg4');
 L5=textread('list_crisSet5_hrg4');
 L6=textread('list_crisSet6_hrg4');
 L7=textread('list_crisSet7_hrg4');

figure(2);clf;plot(L1(:,1),L1(:,2),'bo');grid on;axis([0 730 600 1100]);hold on;
 plot(L2(:,1),L2(:,2),'ro');                                                     
figure(4);clf;plot(L3(:,1),L3(:,2),'bo');grid on;
figure(5);clf;plot(L4(:,1),L4(:,2),'bo');grid on; axis([1590 2240 2150 2560]);hold on;
 plot(L5(:,1),L5(:,2),'ro')
 plot(L6(:,1),L6(:,2),'go')
 plot(L7(:,1),L7(:,2),'mo')

%}   
   
