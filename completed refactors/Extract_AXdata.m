%% Extracts data from AX titrations
% Data is saved in one structure per file
% Created 2020/09/08 MLP based on old extract-data codes
% Updated 2020/11/22 MLP: I have edited the Labview code, so that it
% applies volume correction per every 5 mL.  I also edited it so that on
% each line, the volume listed is what was added, waited 15 s, then pH on
% the corresponding line recorded (i.e., volume and pH value correspond to
% eachother). 
function [w0,S,emf0,t0,CNaOH,CHCl,emf1,t1,w1,emf2,t2,w2,NaOH_ID,HCl_ID,Bur1T,Bur2T,Air1T,Air2T] = Extract_AXdata(file)
% % Open file and extract data
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
w0 = str2double(sample{1})/1000; 
S = str2double(sample{2});
emf0 = str2double(sample{3});
t0 = str2double(sample{4});
CNaOH = str2double(sample{5});
CHCl = str2double(sample{6});
if exist(sample{7}) ~= 0
    NaOH_ID = sample{7};
else
    NaOH_ID = 'nan'; 
end
if exist(sample{8}) ~= 0
    HCl_ID = sample{8};
else
    HCl_ID = 'nan'; 
end

% % See if there is BWD titration data present
Index = find(contains(Lines,'BWD'));
if isempty(Index)
    Index = length(Lines);
    A = 0;
else
    A = 1;
end
% % Extract FWD titration data
clear i;
for i = 2:Index-1
    R         = strsplit(Lines{i,1},',');
    emf1(i-1) = str2double(R{2});
    t1(i-1)   = str2double(R{3});
    V1uncorr(i-1)= str2double(R{6});
    if contains(file, 'sim') == 0
        Bur1T(i-1)= str2double(R{7});
        Air1T(i-1) = str2double(R{9});
    elseif contains(file, 'sim') == 1
        Bur1T(i-1)= 20;
        Air1T(i-1) = 20;
    end
end 
    if contains(file, 'sim') == 0
        V1 = 1.002781*V1uncorr - 0.000200*V1uncorr.^2; %dosimat 12 correction function, updated by MLP 2021/03/10
        w1 = V1.*(1.02888 - 1.069e-4*Bur1T - 4.10e-6*Bur1T.^2);% HCl A21 density function
    elseif contains(file, 'sim') == 1
        w1 = V1uncorr; 
    end
w1 = w1/1000; 
clear R
% % BWD data, if it exists
if A == 0
    emf2 = nan;
    t2   = nan;
    w2   = nan;
else
    for i = Index+1:length(Lines)
    R         = strsplit(Lines{i,1},',');
    emf2(i-Index) = str2double(R{2});
    t2(i-Index)   = str2double(R{3});
    V2uncorr(i-Index)   = str2double(R{6});
    if contains(file, 'sim') == 0
        Bur2T(i-Index)= str2double(R{8});
        Air2T(i-Index) = str2double(R{9});
    elseif contains(file, 'sim') == 1
        Bur2T(i-Index)= 20;
        Air2T(i-Index) = 20;
    end

    end
end
    if contains(file, 'sim') == 0
        V2 = 1.006434*V2uncorr - 0.000267*V2uncorr.^2 ; %dosimat 12 correction function
        rho_NaOH = 1.03121 - 1.172e-4*Bur2T - 4e-6*Bur2T.^2; 
        w2 = V2.*rho_NaOH;% NaOH #J density function 
    elseif contains(file, 'sim') == 1
        w2 = V2uncorr; 
    end
    w2 = w2/1000;
end