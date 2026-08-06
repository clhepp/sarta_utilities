error('this is wrong, called too late in the game ... see interp_aircraft_plevs')

plevs = prof.plevs(:,length(prof.stemp));
plays = plevs2plays(plevs); plays = plays(1:end-1);
palts = prof.palts(:,length(prof.stemp));
dz = abs(diff(palts));

if plevs(1) < plevs(2)
  %% need to flip because kCARTA has TOA at layer 100
  plevs = flipud(plevs);
  plays = flipud(plays);
  palts = flipud(palts);
  dz    = flipud(dz);
end

new_end_palts = max(palts) + max(dz)/1;
new_end_palts = max(palts) + max(dz)/10000;
new_end_plevs = exp(interp1(palts,log(plevs),new_end_palts,[],'extrap'));
new_end_plays = exp(interp1(palts(1:100),log(plays),new_end_palts,[],'extrap'));     
printarray([max(palts)/1000 min(plevs); palts(101)/1000 plevs(101); new_end_palts/1000 new_end_plevs],'in km and mb [max(palts) min(plevs); palts(1) plevs(1); new_end_palts new_end_plevs]')

clear wah
wah.fkc = Y1.fairs;
for iY = 1 : length(mfiles)
  str = ['wah.kc' num2str(iY) ' = squeeze(Y' num2str(iY) '.rairs_all(1,:,:));'];
  eval(str)
end
%wah.kc1  = squeeze(Y1.rairs_all(1,:,:));
%wah.kc2  = squeeze(Y2.rairs_all(1,:,:));
%wah.kc3  = squeeze(Y3.rairs_all(1,:,:));
%wah.kc4  = squeeze(Y4.rairs_all(1,:,:));
%wah.kc5  = squeeze(Y5.rairs_all(1,:,:));

plot(wah.fkc,wah.kc1(:,99:101))
iX = find(wah.fkc >= 660,1); plot(wah.kc1(iX,:),1:101,'+-'); axis([0 1 91 101])
plot(wah.fkc,wah.kc1(:,99)-wah.kc1(:,100),'.-')
figure(5); plot(wah.fkc,wah.kc1(:,101))

if contains(csens,'AIRS')
  for iY = 1 : length(mfiles)
    str = ['Y = Y' num2str(iY) '.rairs_all;'];    
    eval(str);
    Yorig = Y;
    [mm,nn,oo] = size(Y);

    chanschanged(1,nn) = 0;
    for aa = 1 : mm    %% loop over 14 angles
      Yjunk = squeeze(Y(aa,:,:));	   
      for cc = 1 : nn  %% loop over 2378 chans or so
        l2s   = Yjunk(cc,:);
        if l2s(2) < 0.984375  %% there is a big jump to 1.0 transmittance   1-1/64
          chanschanged(cc) = 1;
	  l2s100(cc) = l2s(100);
          l2s(101) = interp1(log(plays(1:100)),l2s(1:100),log(new_end_plays),[],'extrap');
          Yjunk(cc,:) = min(l2s,1);
        end  %% if l2s(2) < 0.984375
      end    %% for cc = 1 : nn
      str = ['Y' num2str(iY) '.rairs_all(aa,:,:) = Yjunk;'];
      eval(str);
    end      %% for aa = 1 : mm

    if iPlot > 0
      figure(1); clf; pcolor(Yorig(:,:,100)); colorbar; colormap jet; shading interp; title('Yorig(:,:,100)')
      figure(2); clf; pcolor(Yorig(:,:,101)); colorbar; colormap jet; shading interp; title('Yorig(:,:,101)')
      figure(3); clf; pcolor(Y1.rairs_all(:,:,100)); colorbar; colormap jet; shading interp; title('YN(:,:,100)')
      figure(4); clf; pcolor(Y1.rairs_all(:,:,101)); colorbar; colormap jet; shading interp; title('YN(:,:,101)')

      % wah.kcN1  = squeeze(Y1.rairs_all(1,:,:));
      % wah.kcN2  = squeeze(Y2.rairs_all(1,:,:));
      % wah.kcN3  = squeeze(Y3.rairs_all(1,:,:));
      % wah.kcN4  = squeeze(Y4.rairs_all(1,:,:));
      % wah.kcN5  = squeeze(Y5.rairs_all(1,:,:));      
      for iY = 1 : length(mfiles)
        str = ['wah.kcN' num2str(iY) ' = squeeze(Y' num2str(iY) '.rairs_all(1,:,:));'];
        eval(str)
      end

      figure(5); plot(wah.fkc,wah.kc1(:,101),'b.-',wah.fkc,wah.kcN1(:,101)); title('level 101 should have changed')
      figure(6); plot(wah.fkc,wah.kc1(:,100),'b.-',wah.fkc,wah.kcN1(:,100)); title('level 100 should be same')
  
      figure(7); hah = find(chanschanged == 1);
      iX = hah(1);              
      iX = find(wah.fkc<= 780); iX = find(l2s100(iX) == min(l2s100(iX)),1);
      iX = find(wah.fkc<= 780); iX = find(l2s100(iX) < 0.2,1);
      if length(mfiles) >= 4
        plot(wah.kc1(iX,:),1:101,'b+-',wah.kcN1(iX,:),1:101,'b--',...
	     wah.kc2(iX,:),1:101,'gs-',wah.kcN2(iX,:),1:101,'g--',...
	     wah.kc3(iX,:),1:101,'r.-',wah.kcN3(iX,:),1:101,'r--',...
	     wah.kc4(iX,:),1:101,'k+-',wah.kcN4(iX,:),1:101,'k--'); axis([0 1 91 101])
      elseif mength(mfiles) < 4
        plot(wah.kc1(iX,:),1:101,'b+-',wah.kcN1(iX,:),1:101,'b--',...
	     wah.kc2(iX,:),1:101,'gs-',wah.kcN2(iX,:),1:101,'g--'); axis([0 1 91 101])

      end	    
    end
    %error('ljkglsjlkjg')
    
  end        %% for iY = 1 : length(mfiles)
else
  error('can only handle AIRS like aircraft right now')
end
clear Y Yjunk Yorig
