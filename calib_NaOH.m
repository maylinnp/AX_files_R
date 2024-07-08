%% Gran Function
% To estimate a) titrant concentration, b) E0, c) protolytic impurities
% (V0+V)*10^-pH.
% As of 2020/03/23, it only calculated c(NaOH). I am working to make this a
% function, so that I can run whatever data I have, it will determine how
% many titrations were performed, and it will calculate an average c(NaOH)
% based on titr 2-end (excluding titration 1). 
% 2020/08/06 I have updated my titration program to include batch
% information, so only NaCl files from #15 and forward fits this algorithm
% now. All NaCl titration files goes to a separate folder. 

calib_NaOH <- function(calibration_files){
  if pwd ~'/Users/May-Linn/Dropbox (Personal)/Matlab/AXdata/NaOH calib';
      cd '/Users/May-Linn/Dropbox (Personal)/Matlab/AXdata/NaOH calib'
  end
  Data = dir(calibration_files);
  titrations = length(Data);
  % % Titration 1---------------------------------
  i = 1;
  file = Data(1).name;
  % [w0,~,w,emf,wacid(i),k,batch]= Extract_data(file);
  [w0,I,~,~,~,CHCl,wacid,emf,t,w,NaOH_ID] = Extract_NaOHdata(file);
  [k,KW,K1,K2,~,~,~,~,~,~,~,~,~,~] = EqConstants (I,mean(t),1,1);
  GF          = (w0+wacid(i)+w).*exp(emf./k);
  GFcurve     = GF(1:length(find(GF > 100)));
  
  figure(1);
  plot(w(1:length(GFcurve)),GFcurve,'o'); hold on; legend('-DynamicLegend')
  F1          = GFcurve;
  ft          = fittype('poly1');
  [fitresult,G(i)]= fit(w(1:length(GFcurve))',F1',ft);
  p1          = fitresult.p1;
  p2          = fitresult.p2;
  weq         = -p1/p2;
  cNaOH(i)    = -p1/p2*CHCl*wacid(i); 
  NaOHex      = (w(end)-(-p2/p1))*cNaOH(i); %in moles, extra base past endpoint 
  HClneutr    = NaOHex/CHCl; %wt of HCl to neutralize extra base
  E0est(i)    = mean(emf(1:length(GFcurve)) - k.*log((wacid(i)*CHCl - w(1:length(GFcurve)).*cNaOH(i))./(w0 + w(1:length(GFcurve)))));
  
  % % Titration 2:end---------------------------------
  for i = 2:titrations
      file = Data(i).name;
  [~,~,~,~,~,~,wacid(i),emf,~,w,~] = Extract_NaOHdata(file);
      w0          = w0 + sum(wacid(1:i-1));
      GF          = (w0+wacid(i)+(w)).*exp(emf./k);
  GFcurve         = GF(1:length(find(GF > 100)));
      figure(1);
      plot(w(1:length(GFcurve)),GFcurve,'o'); hold on; legend('-DynamicLegend')
      F1          = GFcurve;
      ft          = fittype('poly1');
  [fitresult,G(i)]= fit(w(1:length(GFcurve))',F1',ft);
      p1          = fitresult.p1;
      p2          = fitresult.p2;
      cNaOH(i)    = -p1/p2*CHCl*(wacid(i)-HClneutr); 
      NaOHex      = (w(end)-(-p2/p1))*cNaOH(i); %moles of "excess" base, past equivalence point 
      HClneutr    = NaOHex/CHCl; %wt of HCl to neutralize the excess base
      E0est(i)    = mean(emf(1:length(GFcurve)) - k.*log((wacid(i)*CHCl - w(1:length(GFcurve)).*cNaOH(i))./(w0 + w(1:length(GFcurve)))));
  
  end
  mean_cNaOH = mean(cNaOH(2:titrations));
  std_cNaOH  = std(cNaOH(2:titrations))/mean(cNaOH(2:titrations))*100;
  figure(1); title(['NaOH batch ', NaOH_ID]);hold off;
  fprintf('\n\n c(NaOH) = %0.5f mol kg-1 +/- %0.2f %%\n\n\n',mean_cNaOH,std_cNaOH)
  end
  
  %% Sub routines 
  function  [w0,pH,w,emf,wacid,k,batch] = Extract_data(file) %
  % % Extracts data from AX titration csv files, first line has sample info,
  % then vectors from fwd titration, folowed by vectors from bwd titration.
  
  fid = fopen(file);
  % % Create cell array
  i = 1;
  Lines = {};
  while feof(fid) == 0
      Lines{i,1} = fgetl(fid);
      i = i+1;
  end
  fclose(fid);
  % % Extract sample info from first line
  firstLine = Lines{1,1}; 
  sample = strsplit(firstLine,',');
  batch = sample{7};
  w0 = str2double(sample{1})/1000; 
  S = str2double(sample{2});
  P      = strsplit(Lines{2,1},',');
  wacid = str2double(P{4})/1000;
  % % Extract BWD titration data
  clear i;
  for i = 4:length(Lines)
      R      = strsplit(Lines{i,1},',');
      pH(i-3)  = str2double(R{5});
      t(i-3)   = str2double(R{3});
      emf(i-3) = str2double(R{2});
      w(i-3)   = str2double(R{4})/1000;
  end
  k       = mean(8.31451*t/96484.56)

  return(list(cNaOH=cNaOH,E0est=E0est,G=G))
}
