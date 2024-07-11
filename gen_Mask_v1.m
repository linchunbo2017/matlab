clc;
clear;
close all;
mRows = 1032;
nCols = 1792;
%% load light intensity matrix
S1 = load('Normalized_DMD_Intensity.mat');
BP1 = struct2cell(S1);
lightPat = cell2mat(BP1);
lightPat = single(lightPat);
disp('1. load done...');
%% get the scan light intensity matrix
blackLine = zeros(1,1920);
for i=1:24
    lightPat(i,:) = blackLine;
    lightPat(1081-i,:) = blackLine;
end
disp('2. black edge done...');
%% sum intensity, get maximum/minimum and the corresponding index
lightPatIntensitySum = zeros(1032,1792);
% sum in the obliqued direction
lightSum = zeros(1,1792);
for i = 1:mRows   
    shiftNo = floor((i-1)/8);
    lightPatIntensitySum(i,:) = lightPat(i+24,1 + shiftNo :nCols + shiftNo);
end
for i = 1:mRows
    lightSum = lightSum + lightPatIntensitySum(i,:);
end
disp('3. sum done...');
% plot(lightSum);
minVal = min(lightSum);
maxVal = max(lightSum);
[colMin] = find(lightSum==minVal);
[colMax] = find(lightSum==maxVal);
ratio = ( maxVal-minVal)/maxVal;
targetLightIntensity = minVal;
%% sort each column in 
Val = zeros(mRows,nCols);
Ind = zeros(mRows,nCols);
for j=1:nCols
    [Val(:,j), Ind(:,j)] = sort(lightPatIntensitySum(:,j),  'descend');  
end 
disp('4. sort done...');
%% gen the mask matrix
Mask = zeros(mRows,nCols);
MaskLine = zeros(1,nCols);
ValLine = zeros(1,nCols);
IndLine = zeros(1,nCols);
for j=1:nCols
    MaskLine = ones(1,mRows);
    ValLine = Val(:,j)';
    IndLine = Ind(:,j)';
    index = 1;
    for i=1:mRows
        tmpSum = sum(ValLine(1,index:end));
        if (tmpSum > targetLightIntensity)
            MaskLine(IndLine(index)) = 0;
        else
            break;
        end
        index = index + 1;
    end
    Mask(:,j) = MaskLine';
    disp(strcat('5. finish the-',num2str(j),'th column mask '));   
end
 %% 
LightIntensityNew = lightPatIntensitySum .* Mask;
lightSumNew = zeros(1,1792);
for i = 1:mRows   
    lightSumNew = lightSumNew + LightIntensityNew(i,:);
end
minValNew = min(lightSumNew);
maxValNew = max(lightSumNew);
ratioNew = ( maxValNew - minValNew)/maxValNew;
%%
x = 1:1:1792;
figure
plot(x,lightSum,'r');
hold on 
plot(x,lightSumNew,'b');
%%
MaskInDMD = zeros(1080,1920);
for i = 1:mRows
    shiftNo = floor((i-1)/8);
    MaskInDMD(i+24,1 + shiftNo:nCols + shiftNo) = Mask(i,:);
end
save('MaskInDMD.mat','MaskInDMD');
maskImg = mat2gray(MaskInDMD);
%imshow(maskImg);
imwrite(maskImg,'maskImg.bmp')
figure
colormap  Parula;
mesh(MaskInDMD);view(2);xlim([1,1920]);ylim([1,1080]);
%%
disp('done!');
