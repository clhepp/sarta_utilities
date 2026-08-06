clear wah

wah.fkc = Y1.fairs;   %% this one should always exist

fprintf(1,'rejig_AIRS_l2s_for_aircraft.m : iYnum = %2i \n',iYnum)

for iY = iYnum
  %% iYnum is which of the mfiles{1} ...mfiles{5} being worked on
  str = ['wah.kc = squeeze(Y' num2str(iY) '.rairs_all(1,:,:));'];
  eval(str)
end

if iPlot > 0
  figure(8); clf; plot(wah.fkc,wah.kc(:,99:101))
  iX = find(wah.fkc >= 660,1); plot(wah.kc(iX,:),1:101,'+-'); axis([0 1 91 101])
  plot(wah.fkc,wah.kc(:,99)-wah.kc(:,100),'.-')
  figure(5); plot(wah.fkc,wah.kc(:,101))
  pause(0.1)
end

for iY = iYnum
  str = ['Y = Y' num2str(iY) '.rairs_all;'];
  fprintf(1,'  iYnum = %2i str = %s \n',iYnum,str);
  eval(str);
  Yorig = Y;
  [mm,nn,oo] = size(Y);

  chanschanged(1,nn) = 0;
  for aa = 1 : mm    %% loop over 14 angles
    Yjunk = squeeze(Yorig(aa,:,:));	   
    for cc = 1 : nn  %% loop over 2378 chans or so
      l2s   = Yjunk(cc,:);
      if l2s(2) < 0.984375  %% there is a big jump to 1.0 transmittance   1-1/64
        chanschanged(cc) = 1;
        l2s100(cc)      = l2s(100);
        l2s101_orig(cc) = l2s(101);	
        l2s(101) = interp1(log(plays(1:100)),l2s(1:100),log(new_end_plays),[],'extrap');
        Yjunk(cc,:) = min(l2s,1);
        l2s101_new(cc) = l2s(101);		
      end  %% if l2s(2) < 0.984375
    end    %% for cc = 1 : nn
    str = ['Y' num2str(iY) '.rairs_all(aa,:,:) = Yjunk;'];
    eval(str);    
    %fprintf(1,'  lopp iYnum = %2i aa = %2i : %s \n',iYnum,aa,str);    
  end      %% for aa = 1 : mm

  if iPlot > 0
    str = ['Ynew.rairs_all = Y' num2str(iY) '.rairs_all;'];
    eval(str)
    figure(1); clf; pcolor(Yorig(:,:,100)); colorbar; colormap jet; shading interp; title('Yorig(:,:,100)')
    figure(2); clf; pcolor(Yorig(:,:,101)); colorbar; colormap jet; shading interp; title('Yorig(:,:,101)')
    figure(3); clf; pcolor(Ynew.rairs_all(:,:,100)); colorbar; colormap jet; shading interp; title('YN(:,:,100)')
    figure(4); clf; pcolor(Ynew.rairs_all(:,:,101)); colorbar; colormap jet; shading interp; title('YN(:,:,101)')

    for iY = iYnum
      str = ['wah.kcN = squeeze(Y' num2str(iY) '.rairs_all(1,:,:));'];
      eval(str)
    end

    figure(5); plot(wah.fkc,wah.kc(:,101),'b.-',wah.fkc,wah.kcN(:,101)); title('level 101 should have changed')
    figure(6); plot(wah.fkc,wah.kc(:,100),'b.-',wah.fkc,wah.kcN(:,100)); title('level 100 should be same')

    figure(7);
    hah = find(chanschanged == 1);
    iX = hah(1);              
    iX2 = find(wah.fkc <= 780); iX = find(l2s100(iX2) == min(l2s100(iX2)),1); iX = iX2(iX);
    iX2 = find(wah.fkc <= 780); iX = find(l2s100(iX2) < 0.5,1); iX = iX2(iX);
    iX2 = find(wah.fkc <= 780); iX = find(l2s100(iX2) < 0.2,1); iX = iX2(iX);
    plot(wah.kc(iX,:),1:101,'b+-',wah.kcN(iX,:),1:101,'r')
    axis([0 1 81 101]);
    axis([0 1 95 101]);
    title(['iYnum = ' num2str(iYnum) ' chan = ' num2str(iX) ' freq = ' num2str(wah.fkc(iX))])
    pause(0.1)    
  end
  
end        %% for iY = iYnum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write out the stats re:changes to l2s(101)

hah = find(chanschanged == 1);
if length(hah) > 0
  fprintf(1,'iYnum = %2i number of channels whose l2s(101) was changed = %4i \n',iYnum,length(hah))
  fprintf(1,'  layer 100 : %.4f +/- %.4f \n',mean(l2s100(hah)),std(l2s100(hah)))
  fprintf(1,'  layer 101 : orig : %.4f +/- %.4f \n',mean(l2s101_orig(hah)),std(l2s101_orig(hah)))
  fprintf(1,'  layer 101 : new  : %.4f +/- %.4f \n',mean(l2s101_new(hah)),std(l2s101_new(hah)))
  fprintf(1,'  layer 101 : diff : %.4f +/- %.4f \n',mean(l2s101_orig(hah)-l2s101_new(hah)),std(l2s101_orig(hah)-l2s101_new(hah)))
  str = ['Y = Y' num2str(iY) '.rairs_all;'];
  eval(str)
  dada = Yorig(:,:,101) - Y(:,:,101); dada = squeeze(dada); dada = abs(dada(:,hah)); dada = dada(:); fprintf(1,'  layer 101 deltasum = %.12e \n',sum(dada));
  dada = Yorig(:,:,100) - Y(:,:,100); dada = squeeze(dada); dada = abs(dada(:,hah)); dada = dada(:); fprintf(1,'  layer 100 deltasum = %.12e \n',sum(dada));
  dada = Yorig(:,:,099) - Y(:,:,099); dada = squeeze(dada); dada = abs(dada(:,hah)); dada = dada(:); fprintf(1,'  layer 099 deltasum = %.12e \n',sum(dada));  
else
  fprintf(1,'iYnum = %2i number of channels whose l2s(101) was changed = %4i \n',iYnum,length(hah))
end

clear Y Yjunk Yorig chanschanged l2s100 l2s101_orig l2s101_new
pause(0.1)
