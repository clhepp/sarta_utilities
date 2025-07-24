[head hatt prof patt] = rtpread(srcrtp);

IPROF=52;

fnames=fieldnames(prof)
fwant={'plevs','gas_';'ptemp'}
iif = [];
iif = contains(fnames,'ptemp');
iif = or(iif, contains(fnames,'gas_'));
iif = or(iif, contains(fnames,'plevs'));

savVars = fnames(iif);

junk = zeros(101,length(savVars));
for ig = 1:length(savVars)
 junk(:,ig) = prof.(savVars{ig})(:,IPROF);
end

FH = fopen('atmos_52.txt','w')
for il = 1:101
  fprintf(FH,'%e %f7.3 %e %e %e %e %e %e %e %e %e\n', junk(il,:));
end

fclose(FH)


prof.stemp(IPROF) = 240.1600
prof.emis(:,IPROF)=
    0.9767
    0.9829
    0.9889
    0.9919
    0.9925
    0.9913
    0.9891
    0.9875
    0.9863
    0.9851
    0.9840
    0.9831
    0.9776
    0.9770
    0.9764
    0.9756
    0.9745
    0.9732
    0.9715
prof.efreq(:,IPROF) = 
1.0e+03 *
    0.7692
    0.8000
    0.8333
    0.8696
    0.9091
    0.9524
    1.0000
    1.0417
    1.0870
    1.1364
    1.1905
    1.2500
    2.4390
    2.5000
    2.5641
    2.6316
    2.7027
    2.7778
    2.8571
prof.satzen(IPROF) =   17.9223
prof.solzen(IPROF) =   150

