%% Extracts data from AX titrations
% Data is saved in one structure per file
% Created 2020/09/08 MLP based on old extract-data codes
% Updated 2020/11/22 MLP: I have edited the Labview code, so that it
% applies volume correction per every 5 mL.  I also edited it so that on
% each line, the volume listed is what was added, waited 15 s, then pH on
% the corresponding line recorded (i.e., volume and pH value correspond to
% eachother). 
function [w0,I,emf0,t0,CNaOH,CHCl,wHCl,emf,t,w,NaOH_ID] = Extract_NaOHdata(file)
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
I = str2double(sample{2});
emf0 = str2double(sample{3});
t0 = str2double(sample{4});
CNaOH = str2double(sample{5});
CHCl = str2double(sample{6});
if exist(sample{7}) ~= 0
    NaOH_ID = sample{7};
else
    NaOH_ID = 'nan'; 
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
    R         = strsplit(Lines{2,1},',');
    wHCl = str2double(R{6})/1000;

clear R 
% % BWD data, if it exists
if A == 0
    emf = nan;
    t   = nan;
    w   = nan;
else
    for i = Index+1:length(Lines)
    R         = strsplit(Lines{i,1},',');
    emf(i-Index) = str2double(R{2});
    t(i-Index)   = str2double(R{3});
    V2uncorr(i-Index)   = str2double(R{6});
    Bur2T(i-Index)= str2double(R{8});
    end
end
    V2 = 1.006434*V2uncorr - 0.000267*V2uncorr.^2; %dosimat 12 correction function
    rho_NaOH = 1.03121 - 1.172e-4*Bur2T - 4e-6*Bur2T.^2; 
    w = V2.*rho_NaOH;% NaOH #J density function 
    w = w/1000;
end