


fin='/home/chepplew/gitLib/ftc_dev/run/jout_set1_fowp';
FH=fopen(fin,'r')

clear ichan npred rms
k = 1;
n = 1;
while ~feof(FH)
  cline = fgetl(FH);

  if ( ~isempty(strfind(cline(1:2),'!')) )
      continue
  end
  if ( isempty(strfind(cline(1:2),'!')) && ~isempty(strfind(cline, 'ichan')) )
      % record channel number and reset n
      ichan(k) = str2num(cline(8:end));
      k = k + 1;
      n = 1;
      continue
  end
  if ( isempty(strfind(cline(1:2),'!')) && ~isempty(strfind(cline, 'npred')) )
     % record npred and rms
     junk = strsplit(cline, ' ');
     npred(k-1,n) = str2num(junk{4});
     rms(k-1,n)   = str2num(junk{7});
     n = n + 1;
     continue;
  end
  %if ( k == 3 ) break;  end
end

