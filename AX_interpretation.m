%% Simple AX interpretation code
% 2021/01/29 Based on AX_interp_v21, I've made a function to go with
% simlated or real data. It reads csv files and extracts data using
% Extract_AXdata.m. Will read CHCl and CNaOH from FileName, unless other
% concentrations are explicitly stated. 

% 2021/02/01 Edited the interp1 and interp2 to utilize incoming [H+] on the
% total scale (if ST is 0 then total scale is the same as free scale). Have
% not added this to the other part of the processing code yet, including
% the plotting part. -MLP

% 2021/02/02 This code works for Gran and nonlin calc of E0 and AT for FWD
% and BWD data (low pH range), although E0 comes out a touch high from BWD
% data. Like there is acid missing (bc it was low when there was extra acid
% from CT). -MLP
% Needed calculations: if there is silicate in the sample, and I want to
% calculate KW, how can I correct for that? 

% 2021/04/19 Major problem solved! Had a typo for narrowing chopping data
% for BWD processing: had used E01est instead of the correct E02est, which
% is why FWD processing pH range affected the quality of the BWD data
% processing. Now that this is rectified, I can perfectly represent the
% data at low pH and the rest of the curve too. I think I might be ready
% for some real live titrations!? -MLP

% NOTE TO SELF: Should set this up to interpret data even without FWD
% titration data, like the NaOH files where one big slew of acid is added.
% This can cut down on some analysis time if I want to

function [Results,ResultsText] = AX_interpretation(FileName)
global H m m02 m0 FT KF CHCl CNaOH KW mHCl ST KS S E02est k AT BT KB
% Extract data from FileName
if pwd ~'/Users/May-Linn/Dropbox (Personal)/Matlab/AXdata';
    cd '/Users/May-Linn/Dropbox (Personal)/Matlab/AXdata'
end

% Check sample type 
    if contains(FileName,'NaCl') == 1
        Sol = 1;
    elseif contains(FileName,'KCl') == 1
        Sol = 3;
    elseif contains(FileName,'Cl') == 0
        Sol = 2; 
    end

% Air buoyancy correction if it is a real file, bc _sim files are practically in mass    
    if contains(FileName, '.csv') == 1
        [w0,S,emf0,t0,CNaOH,CHCl,emf1,t1,w1,emf2,t2,w2,NaOH_ID,HCl_ID,Bur1T,Bur2T,~,~] = Extract_AXdata(FileName);
        rho_HCl  = 1.02888 - 1.069e-4*Bur1T - 4.10e-6*Bur1T.^2;
%         rho_NaOH = 1.03121 - 1.172e-4*Bur2T - 4e-6*Bur2T.^2; %#J-2
        rho_NaOH = 1.02727 - 1.088E-04*Bur2T -4.00E-06*Bur2T.^2; %batch #K-2, average 
        rho_Sol  = 1.027;
        m1 = w1.*((1-(0.0012013/8))./(1-(0.0012013./rho_HCl)));
        m2 = w2.*((1-(0.0012013/8))./(1-(0.0012013./rho_NaOH)));
        m0 = w0*((1-(0.0012013/8))/(1-(0.0012013/rho_Sol)));
    elseif contains(FileName, '.mat') == 1
        load(FileName);
    end 
    m02     = m0 + m1(end); 
    BTratio = 1;
    CT_NaOH = 0e-6;
    SiT_NaOH = 0e-6; 
    if Sol == 1 || Sol == 3
       CT_sample = 2.5e-6; 
       SiT = 1e-6;
       PT = 0e-6;
    else
       CT_sample = 4e-6; 
       SiT = 5e-6;
       PT = 0.5e-6;
    end
    
    
    [k,KW,K1,K2,ST,FT,BT,KS,KF,KB,KSi,KP1,KP2,KP3] = EqConstants(S,t2,Sol,BTratio); 
    options = optimset('TolX',1e-32,'TolFun',1e-32,'Display','off','Algorithm','levenberg-marquardt');
    pH1est  = -log10(exp((emf1 - 0.41)./k)); %this E0 guess is important 
    pH2est  = -log10(exp((emf2 - 0.41)./k));
    Idx1    = find (pH1est > 3 & pH1est < 3.5); %pH 3.5–3
    Idx3    = find (pH2est >  3 & pH2est < 3.5); %pH 3.5–3
    Idx2    = find (pH2est < 10.5 & pH2est > 9); %pH 9.5-11 or so

%% Here I want to add FWD interpretation, if such data exist in FileName
if length(m1) > 1
    [p11FWD,p21FWD,weq1FWD] = Gran1(m0,m1,emf1,Idx1);
    E01est      = (emf1(Idx1) - k.*log(CHCl.*(m1(Idx1) - weq1FWD)./(m0 + m1(Idx1))));
    H           = exp((emf1 - mean(E01est))./k); 
    m           = m1; 
    AT1est      = -p21FWD/p11FWD*CHCl/m0; 
    % Non-linear solver
    x0(1)       = 1;     %f
    x0(2)       = AT1est; %AT
    [x]         = fsolve(@interp1, x0, options);
    r1          = 1000000*interp1(x); 
    sos_1       = r1*r1';
    E01         = mean(E01est) - k*log(x(1));
    H           = exp((emf1 - E01)./k); 
    AT_FWD      = x(2); 
else
    E01 = nan; 
    AT_FWD = nan; 
end

%% Data processing of BWD data, including E0 and, if -Cl solution, KW 
% Low-pH side of data, estimate E0 and AT using Gran
    [p11BWD,p21BWD,weq1BWD] = Gran1(m02,m2,emf2,Idx3);
    E02est  = (emf2(Idx3)-k.*log(CNaOH.*(weq1BWD - m2(Idx3))./(m02+m2(Idx3))));
    E0      = mean(E02est); 
    AT2est  = CHCl*m1(end)/m0 + (p21BWD/p11BWD*CNaOH)/m0; 
% Low-pH side of data, estimate E0 and AT using non-linear 
    x0(1)       = 1;
    x0(2)       = AT2est; 
    if length(m1) > 1
        AT          = AT_FWD;
    end
    H           = exp((emf2(Idx3)-mean(E0))/k); 
    m           = m2(Idx3);
    mHCl        = m1(end);
    [x]         = fsolve(@interp2,x0, options);
    r2          = 1000000*interp2(x); 
    sos_2       = r2*r2';
    E02         = (mean(E02est) - k*log(x(1)));
    if length(m1) > 1
        H1T         = exp((emf1-E01)./k);
    end
    AT_BWD      = x(2);

if Sol == 1 || Sol == 3 
 % Alkaline data, calculate KW (for synthetic solutions only)
    [p12BWD,p22BWD,weq2BWD,KWest,KWstd] = Gran2(m02,m2,emf2,Idx2);
    Xest       = (weq2BWD-weq1BWD)*CHCl/m0;
    KW = KWest; 
end 


%% AX curve expression 

    H2T    = exp((emf2-E02)./k) ;
    pHT    = -log10(H2T);
    BOH4   = m0*(BT)./(1+(H2T./KB));
    if contains(FileName, 'sim') == 0
        pKW    = -log10(KW) ; %** This correction comes from Dickson and Riley, comment on Hansson's determination
        KW     = 10^-pKW;    
    end
    Si     = (m0*SiT + m2.*SiT_NaOH)./(1+(H2T./KSi));
    HF     = m0*FT./(1+KF./H2T);
    Z      = 1 + ST/KS; 
    HSO4   = m0*ST./(1+(Z*KS)./(H2T)); 
    AC     = (m2.*CT_NaOH + m02*CT_sample).*((K1.*H2T+2*K1*K2)./(H2T.^2 + K1.*H2T + K1*K2)); 
    AP     = m0*PT.*((KP1.*KP2.*H2T + 2.*KP1.*KP2.*KP3 - H2T.^3)./(H2T.^3 + KP1.*H2T.^2 + KP1.*KP2.*H2T + KP1.*KP2.*KP3));
    
    logKMgOH2 = 2.8; KMgOH2 = 10^(-logKMgOH2); 

    CMg = 50e-3; %nozaki table value
    logKMgOH  = 350; KMgOH  = 350;%10^-logKMgOH;
    % K = [Mg][OH]/[MgOH], [MgOH] = K/([Mg][OH])
    OH = Z.*(KW)./H2T; MgOH2 = (CMg.*OH.^2)./KMgOH2;
    MgOH = 0;%(CMg.*OH).*KMgOH;
    
    AXcurve= (m0*AT_BWD - CHCl*mHCl + m2.*CNaOH + HF + HSO4 - BOH4 - Si - AC - AP + ...
        (m0 + mHCl + m2).*((H2T./Z) - (Z.*(KW)./H2T)) + MgOH + 2.*MgOH2)./(m0);
    
    figure(1)
    hold on
    plot(pHT,AXcurve*1e6,'s-')
    ylabel('\Delta{\itA}_{X} (\mumol kg^{-1})')
    xlabel('pH_{T}')
% 
% 
%     xlim([2 9])
%     ylim([-10 inf])
    
%% Fitting the data to stuff
H0  = (exp((emf0-E02)./k));



% what will this pH be after CT is moved to basically 0 umol/kg? a tad
% higher? 
% look for this pH in the AX curve. What if not spot on? How do I
% interpolate? 
Idx4 = find(H2T > H0);
dHfromH0 = (H2T(Idx4(end)) - H0)/(H2T(Idx4(end)) + H2T(Idx4(end)+1)); % Distance to H0 from the closest lower [H+]
AXtotal = AXcurve(Idx4(end))*(1+dHfromH0);

ResultsText = sprintf('For sample %s: \rE0 estimated from FWD data is %.6f V. \rE0 estimated from BWD data is %.6f V. \rAT estimated from FWD is %.2f umol/kg. \rAT estimated from BWD is %.2f umol/kg. \rAXTotal is %.2f umol/kg\r\n',FileName,E01,E02,AT_FWD*1e6,AT_BWD*1e6, AXtotal*1e6);

% % when calculating the euqilibrium pH, the points are probably far away
% bc pH changes quickly here. Therefore, I must change both points by + and
% -, respectively, the 1 stdev in emf and see how big of a change that
% produces. That will be my error in estimating pH. weq probably doesn't
% matter there as it is a small portion of the overall, and each weight
% step is the same. 

Results.AT_BWD = AT_BWD;
Results.E0_BWD = E02;
if length(m1) > 1
    Results.AT_FWD = AT_FWD;
    Results.E0_FWD = E01;
    Results.m1     = m1;
    Results.H1T    = H1T;
else
    Results.AT_FWD = AT_FWD;
    Results.E0_FWD = E01;  
    Results.m1     = nan;
    Results.H1T    = nan;
end
Results.AXtotal= AXtotal;
Results.m2     = m2
Results.H2T    = H2T; 
Results.AXcurve= AXcurve; 
Results.emf2   = emf2; 

%% Identifying peaks above uncertainty 
% Basically, until pH 8 that is below 2.5 umol/kge

% Re-calculate "AT pK equilibrium value" under current T/S conditions
% At pH values below 8.something where uncertainties really rise:
% Look for slopes that span more than 2.5 umol/kg from bottom to top, any
% slopes that span less than 2.5 umol/kg in AXcurve are not to be counted 
% (this is all before I include silicate in the analysis, and uncertainty in 
% that value). 

% % % 1. Look for change between two points that is larger than 2.5 umol/kg in AX.
% % % 2. If two points has a slope larger, add these two to a new string of indexes.
% % % 3. Keep checking the slope between the following points, once the slope is 
% % % below 2.5, remove that last point, use that set of indexes to produce a slope
% % % from which KX can be calculated, and where the top point approximates [X] 
% % % (but not XT). So the last index + 1 is used to estimate X, and with KX we 
% % % can calculate XT too. 
% % % 4. It makes sense to limit the span of this "looking for slope" to within 
% % % 3 pH units, or within 3 orders of magnitude, because beyond that I don't 
% % % have great enough sensitivity to measure changes in X. If the slope continues,
% % % it is likely another X. If two Xs are overlapping that is a different issue. 
% % % 
% % Nonlinear fit of data, no constraints 
% t = fittype('XT/(1+(H2T/KX)) + X1T/(1+(H2T/KX1))','coefficients', {'XT','KX','X1T','KX1'},'indep','H2T');
% A = find(AX > 0.99*max(AX));
% B = find(pHT > 3);
% s = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0,10^-14,0,10^-14],...
%                'Upper',[max(AX), 10^-4.5, max(AX), 10^-4.5],...
%                'Startpoint',[AX(A(end)), 10^-8.5, AX(A(end))/2, 10^-4.5]);
% fitAX = fit(H2T(B:A(end))',AX(B:A(end))',t,s)

% % Nonlinear optimization of data, with constraints
% X1T = optimvar('X1T');
% KX1 = optimvar('KX1');
% % X2T = optimvar('X2T');
% % KX2 = optimvar('KX2');
% 
% prob = optimproblem;
% prob.Objective = X1T./(1+(H2T./KX1)) ;
% prob.Constraints.cons1 = X1T + X2t <= max(AX);
% prob.Constraints.cons2 = X1T >= 0;
% % prob.Constraints.cons3 = X2T >= 0;
% x0.X1T = max(AX);
% x0.KX1 = 10^-9.5;
% % x0.X2T = 0;
% % x0.KX2 = 10^-4.5;
% sol = solve(prob,x0)
% 
% figure(1)
% hold on
% plot(pHT(B:A(end))',fitAX(H2T(B:A(end)))*1e6,'^:')
end 
%% Interpretation function
function r1 = interp1(x) 
global m0 CHCl m H KW KF KS ST FT 
    f    = x(1); 
    AT   = x(2);
    Z    = 1 + ST/KS; 
    HSO4 = m0.*ST./(1+(Z*KS)./(f*H)); 
    HF   = m0.*FT./(1+KF./(f*H));
    r1   = m0*AT + HSO4 + HF - m.*CHCl + (m0 + m).*((f*H./Z) - (Z.*KW./(f*H)));
end

function r2 = interp2(x) 
global m H KW KF FT CHCl mHCl m0 CNaOH KS ST
    f    = x(1); 
    AT   = x(2);
    Z    = 1 + ST/KS; 
    HSO4 = m0.*ST./(1+(Z*KS)./(f*H)); 
    HF   = m0.*FT./(1+KF./(f*H));
    r2   = ((m0*AT + HSO4 + HF - CHCl*mHCl + m.*CNaOH + (m0 + mHCl + m).*((f*(H./Z)) - (Z.*KW./(f*(H))))))./(m0 + mHCl + m);
end

function [p1,p2,weq] = Gran1(m0,m,emf,Idx)
global k
    F1          = (m0+m(Idx)).*exp((emf(Idx))./k);
    ft          = fittype('poly1');
    fitresult   = fit(m(Idx)',F1',ft);
    p1          = fitresult.p1;
    p2          = fitresult.p2;
    weq         = -p2/p1;
end

function [p1,p2,weq,KWest,KWstd] = Gran2(m0,m,emf,Idx)
global k E02est CNaOH
    F2         = (m0+m(Idx)).*exp(((-emf(Idx))./k));
    ft         = fittype('poly1');
    fitresult  = fit(m(Idx)',F2',ft);
    p1         = fitresult.p1;
    p2         = fitresult.p2;
    weq        = -p2/p1;
    KW         = (CNaOH.*(m(Idx) - weq)./(m0 + m(Idx)))./(exp((-emf(Idx) + mean(E02est))./k)); 
    KWest      = mean(KW);
    KWstd      = std(KW); 
end 