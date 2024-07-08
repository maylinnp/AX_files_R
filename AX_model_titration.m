%% AX titration model

% 2021/01/21 This script models AX titrations, and realistic uncertainties
% will be added to various parts of the titration. This includes scatter un
% the dosimat, and carbonate and silicate contamination of the titrants
% (primarily NaOH). The data output will be formatted the same as the
% output files from my Labview program

% 2021/01/26 Working on the titration part of the code (not th equilibrium
% so far). I need to figure out a better naming system for the two weight
% vectors. 

% 2021/01/27 I believe I have the whole code working, last thing is
% figuring out the file naming structure which should be part of the
% function input. Smaller increment size is needed to simulate -Cl solutions! 

% 2021/01/28 There is an issue with the code. Especially at low pH, when I
% interpret my data (BWD), there is an obvious scatter that diminishes at
% higher pH values. Probably because the change is very small for each
% increment? I don't know, anyways, when I interpret the data the scatter
% becomes clear. I think it has to do with what solver I use and how I set
% it up, like ther eare just so many solutions and guess/number of
% iterations matters. Maybe I can make it a scale problem? Not sure. Brain
% tired. Will look at it tomorrow. There is currently no other noise added
% into the system, so it must be coming from the interpretation. -MLP

% 2021/01/29 The scatter at low pH was NOT caused by the solver, it was
% caused by me saving the variables in .csv file and specifying only 6
% decimal points. Once I upped this to 15, the scatter disappeared!
% (Thanks Mike!). 

% 2021/02/01 I added for H+ in the ChemEq to be on the total scale, meaning
% that the appropriate acid-bases are corrected to free-scale if their
% constants require it, but most are not. This by adding a factor, Z,
% according to Dickson et al. 2003. This should work for NaCl and KCl
% solutions as well, because their ST is zero and Z will therefor be 1.
% -MLP

% 2021/02/02 Finally fixed the issue that was interpretation code not
% giving me the AT and E0 I put in – in my titration code I had reduced CT
% from initial to "removed" value one line too late, so the first titration
% point still had 2000 umol/kg CT left in it. I could see this even on the
% gran plot, so I should have understood that the issue was related to the
% first 3.5 titration point. Oh well, de-gassed CT was moved to the
% appropriate place and I can keep working on the BWD parts of the two
% codes. -MLP

% 2021/04/23 I want to update my code so I don't input AT, but rather sums
% up AT and excludes AX? Or, what do I really need here. For now perhaps
% make the code so it's quick and easy for me to use. Worry about adding
% user friendly stuff later. I do want it as a function though, and I'd
% like to do something similar with the interpretation code. 

% 2021/5/09 Started adding uncertainties in a random way (I hope) so I can
% run many in a Monte Carlo type simulation, and get an average and
% standard deviation for the LOD and "0 AX signal". -MLP

function fileNameGen = AX_model_titration(sample_type,S_or_I,samplenumber)%,Sample,System,AX,Uncert)

% arguments
%     sample_type = 'SW'
%     S_or_I = 35
%     samplenumber = 1
%     Sample.AT = 2200e-6 
%     Sample.CT = 2000e-6
%     System.CNaOH = 0.05
%     System.CHCl = 0.1
%     System.Temp = 20
%     AX.X1T = 0e-6
%     AX.pKX1 = 5
%     AX.X2T = 0e-6
%     AX.pKX2 = 8
%     Uncert.emf_noise = 0e-6
%     Uncert.w_noise = 0.00000e-3
%     Uncert.T_noise = 0.0
%     Uncert.E0_drift = 0.000
%     Uncert.BTrat = 1
% end
global d k KW K1 K2 ST FT BT KS KF KB KSi KP1 KP2 KP3 KX1 KX2 X1T X2T SiT PT CT AT 
% % Optional function inputs, will be set to standard value defined below
% if not defined during function call 
    

    CT = 2000e-6;
    CHCl = 0.1 ;
    CNaOH = 0.05 ;
    t = 20;
    X1T = 0E-6;
    KX1 = 10^-5;
    X2T = 0E-6;
    KX2 = 10^-7;
    Sample.AT = 2200e-6 + ( X1T + X2T);
    emf_noise = 15e-6; %where does this number come from 
    w_noise = 0.0007e-3; 
    T_noise = 0.0;
    E0_drift = 0;
    BTrat = 1;

%% Set sample conditions 
T0      = 273.15 + t;
if strcmp(sample_type,'NaCl') == 1
    Sol = 1;
elseif strcmp(sample_type,'KCl') == 1
    Sol = 3;
elseif strcmp(sample_type,'SW') == 1
    Sol = 2;
end

[k,KW,K1,K2,ST,FT,BT,KS,KF,KB,KSi,KP1,KP2,KP3] = EqConstants(S_or_I,t,Sol,BTrat); 

% add uncertainty according to sources presented in DOE/Dickson (2007)
pK1 = -log10(K1) - 0.007; K1 = 10^-pK1;
pK2 = -log10(K2) - 0.01; K2 = 10^-pK2; 
pKW = -log10(KW) - 0.007; KW = 10^-pKW;
pKF = -log10(KF) + 0.02; KF = 10^-pKF;
% pKS = -log10(KS) - 0.1; KS = 10^-pKS;
pKB = -log10(KB) - 0.004; KB = 10^-pKB;
pKSi = -log10(KSi) - 0.02; KSi = 10^-pKSi;
pKP1 = -log10(KP1) + 0.09; KP1 = 10^-pKP1;
pKP2 = -log10(KP2) + 0.03; KP2 = 10^-pKP2;
pKP3 = -log10(KP3) - 0.2; KP3 = 10^-pKP3;

wa      = 0.05e-3;    % increment HCl
wb      = 0.025e-3;    % increment NaOH
w0real  = 115e-3;     % sample weight
w0      = w0real ;%+ 0.15e-6*randn;
E0      = 0.4;        % true E0 of the electrode 
pHend   = 10;        % endpoint for BWD titration
PT      = 1e-6;

% Sample contamination 
CT_degas    = 4e-6;        % CT left in sample after degassing                   
CT_NaOH     = 0e-6;        % CT in NaOH titrant
SiT_sample  = 15e-6;        % these will be added in the same way as CT contamination
SiT_NaOH    = 0e-6; 
SiT_HCl     = 0e-6; 
%% Calculate initial conditions and FWD titrate 
options = optimset('TolX',1e-32,'TolFun',1e-32,'Display','off','Algorithm','levenberg-marquardt');
H_guess = 10^-8; % could change this for CO2sys 
d       = 1;
AT      = Sample.AT;
SiT     = SiT_sample;
[H0]    = fzero(@ChemEq,H_guess,options); % I might want to change this to fsolve at some point 
emf0    = log(H0)*k + E0;

% Add acid to reach pH 3.5 to de-gas 
H_guess = 10^-3.5;
i       = 1;
CT      = CT_degas; 
wHCl(i) = (-H_guess*w0 - Sample.AT*w0) / (H_guess - CHCl); 
wHCl_real(i)= wHCl(i)+w_noise*randn(1,1); 
AT      = (Sample.AT*w0 - wHCl_real(i)*CHCl) / (w0 + wHCl_real(i));
d       = w0 / (w0 + wHCl_real(i));
[H1(i)] = fsolve(@ChemEq,H_guess,options);

% Titrate from pH ~3.5 to ~3
while H1(i) < 10^-3
  i       = i + 1;
  wHCl(i) = wHCl(i-1) + wa; 
  wHCl_real(i)=wHCl_real(i-1) + wa + w_noise*randn(1,1); 
  d       = w0/(w0 + wHCl_real(i)); 
  AT      = (Sample.AT*w0 - CHCl*wHCl_real(i)) / (w0 + wHCl_real(i)); 
  H_guess = H1(i-1);   
  [H1(i)] = fsolve(@ChemEq,H_guess,options); 
end
% Convert data and add noise
emf1= (log(H1).*k + E0) + emf_noise*randn(length(H1),1)';
pH1 = -log10(H1); T1 = ones(1,length(wHCl),1)*t + T_noise*randn(1,length(wHCl)); 
%% BWD titration 
E0      = E0 + E0_drift; % to mimic any drift in this parameter that could follow the 40 minute de-gassing period
i       = 1;
w02     = w0 + wHCl_real(end);
wNaOH(i)= 0; wNaOH_real(i) = wNaOH(i);
emf2(i) = log(H1(end))*k + E0;
H2(i)   = exp((emf2(1)-E0)/k); % this value will change if there has been drift in E0
ATex    = (Sample.AT*w0 - CHCl*wHCl_real(end))/w02;
while H2(i) >= 10^-pHend
  i     = i+1;
  wNaOH_real(i) = wNaOH_real(i-1) + wb + w_noise*randn(1,1);
  wNaOH(i)= wNaOH(i-1) + wb;
  d     = w0 / (w02 + wNaOH_real(i));
  AT    = ((ATex*w02 + CNaOH*wNaOH_real(i)) / (w02 + wNaOH_real(i)));
  CT    = (CT_degas*w02 + CT_NaOH*wNaOH_real(i))/(w02 + wNaOH_real(i));
  SiT   = (SiT_sample*w0 + SiT_HCl*wHCl_real(end) + SiT_NaOH*wNaOH_real(i))/(w02 + wNaOH_real(i));
  H_guess= [10^-(-log10(H2(i-1))+4) H2(i-1)];
  [H2(i)]= fzero(@ChemEq,H_guess,options);
end
emf2 = (log(H2).*k + E0) + emf_noise*randn(length(H2),1)';
pH2  = -log10(H2); T2 = ones(1,length(wNaOH))*t + T_noise*randn(1,length(wNaOH)); 
%% Save file as .csv file 
% Convert output and save as .csv file, output file has same format as LabView file, minus date/time
DateTime = datetime('now','Format','yyyyMMdd');
fileName = sprintf('%s %s-%d_sim.mat', DateTime, sample_type,samplenumber);
fileNameGen= sprintf('%s %s-*_sim.mat', DateTime, sample_type);
% % I'd like to add a uigetdir command here, so where the user can add a
% preferred file path, and if the path has already been set it doesn't have
% to change but I need to find a way they can change it 
if pwd ~'/Users/May-Linn/Dropbox (Personal)/Matlab/AXdata';
    cd '/Users/May-Linn/Dropbox (Personal)/Matlab/AXdata'
end
% 
% fid = fopen(fileName,'w');
% fprintf(fid,'%.15f, %.15f, %.15f, %.15f, %.6f, %.6f, %s, %s \r\n', w0real*1000, S_or_I, emf0, T0, CNaOH, CHCl, 'B','A1-1');
% for i = 1:length(emf1)
%     fprintf(fid,'%d, %.20f, %.15f, %.15f, %.15f, %.15f \r\n', i, emf1(i), T1(i), wHCl(i)*1000, pH1(i), wHCl_real(i)*1000);
% end
% fprintf(fid,'%s \r\n','BWD');
% for i = 1:length(emf2)
%     fprintf(fid,'%d, %.20f, %.15f, %.15f, %.15f, %.15f \r\n', i, emf2(i), T2(i), wNaOH(i)*1000, pH2(i), wNaOH_real(i)*1000);
% end
%     fclose(fid);
m0 = w0*1000; m1 = wHCl*1000; m2 = wNaOH * 1000; 
S = S_or_I; t0 = T0; t1 = T1; t2 = T2; 
save(fileName,'m0', 'S', 'emf0', 't0', 'CNaOH', 'CHCl', 'emf1', 't1', 'm1', 'pH1', 'emf2', 't2', 'm2', 'pH2')
end %of AX_model
%% Additional functions 
% Chemical equilibrium model used where [H+] is calculated by minimizing
% the residual
function r = ChemEq(x)
global d AT KW KX1 X1T KX2 X2T CT K1 K2 KS ST KB BT KF FT KP1 KP2 KP3 KSi SiT PT
    H    = x;
    Z    = 1 + ST/KS;
    OH   = KW/(H/Z);
    HSO4 = d*ST/(1 + (KS*Z)/H); 
    HF   = d*FT/(1 + KF/H);
    BOH4 = d*BT/(1 + H/KB);
    HCO3 = d*CT*K1*H/(H^2 + K1*H + K1*K2); % check the use of d 
    CO3  = K2*HCO3/H;
    X1   = d*X1T/(1 + H/KX1);
    X2   = d*X2T/(1 + H/KX2);
    H2PO4= d*(PT*KP1*H^2)/(H^3 + KP1*H^2 + KP1*KP2*H + KP1*KP2*KP3);
    HPO4 = KP2*H2PO4/H;
    PO4  = KP3*HPO4/H;
    H3PO4= H2PO4*H/KP1;
    P    = (2*PO4 + HPO4 - H3PO4);
    Si   = SiT/(1 + H/KSi);
    
    r    = (H/Z - OH) + AT + HSO4 + HF - HCO3 - 2*CO3 - BOH4 - X1 - X2 - P - Si;
end