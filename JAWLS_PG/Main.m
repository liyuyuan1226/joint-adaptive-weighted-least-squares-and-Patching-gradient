%*********************************************************************************************************************
%===========Author: Yuyuan Li=========================================================================================
%===========paper: <Depth super-resolution using joint adaptive weighted least squares and Patching gradient>=========
%===========base code from: Liu, Wei, et al. "Robust Color Guided Depth Map Restoration."=============================

%===========The paper has been published on ICASSP 2018===============================================================
clear;
close all;

Depth = double(imread('E:\Robust-Color-Guided-Depth-Map-Restoration-master\Input & Results\Input\art_big.png'));
Color = double(imread('E:\Robust-Color-Guided-Depth-Map-Restoration-master\Input & Results\Input\art_color.png'));
% Color = double(imread('E:\smooth\code\result\GFenhance.png'));
DepthNoise = double(imread('E:\Robust-Color-Guided-Depth-Map-Restoration-master\Input & Results\Input\art_depth_3_n.png'));
% DepthNoise=Depth;%downsample
[m, n] = size(Depth);
% DepthNoise = imresize(DepthNoise, [m/8, n/8]);%downsample
DepthNoise = imresize(DepthNoise, [m, n]);
% DepthNoise=uint8(DepthNoise);
% imwrite(DepthNoise,'E:\ICASSIP2018\my paper\patchgradient\DepthNoise.png');
% patch1=DepthNoise(841:1080,141:540);
% patch2=DepthNoise(141:380,840:1240);
% imwrite(patch1,'E:\ICASSIP2018\my paper\patchgradient\patch1_DepthNoise.png');
% imwrite(patch2,'E:\ICASSIP2018\my paper\patchgradient\patch2_DepthNoise.png');

edge_Sobel=edge(DepthNoise,'Sobel',1);

% edge_Sobel1=zeros(m, n);
% edge_Sobel1(edge_Sobel==0)=1;
figure(431);imshow(edge_Sobel);
% imwrite(edge_Sobel1,'E:ICASSIP2018\my paper\patchgradient\Sobel.png');
% patch1=edge_Sobel1(841:1080,141:540);
% patch2=edge_Sobel1(141:380,840:1240);
% imwrite(patch1,'E:\ICASSIP2018\my paper\patchgradient\patch1_SOBEL.png');
% imwrite(patch2,'E:\ICASSIP2018\my paper\patchgradient\patch2_SOBEL.png');

%自加测试8x
% DepthNoise2 = double(imread('.\Images\depth_3_n.png'));
% DepthNoise2 = imresize(DepthNoise2, [m, n]);
% DepthUpdate2 = DepthNoise2;
% RS2 = GetRelSmooth(DepthUpdate2, 1, 1);
% figure(330);imshow(RS2);
% BW2 = ones(m, n);
% BW2(RS2<0.96)=0;
% figure(331);imshow(BW2);
% 
% BW2 = ones(m, n);
% BW2(RS2<0.9)=0;
% figure(332);imshow(BW2);
% 
% BW2 = ones(m, n);
% BW2(RS2<0.8)=0;
% figure(333);imshow(BW2);

%over


alpha = 0.95;

rSmooth = 9;
sigmaSSmooth = 9;

rRelSmooth = 3;%原参数为7
sigmaC = 10;


BWAdp = [10, 4, 7, 10];

IterNumFixedPoint = 16; % 4/8/16/32 for 2X/4X/8X/16X upsampling

DepthUpdate = DepthNoise;
% DepthUpdate=double(imread('E:\result\result_test2\result_24.png'));
t11 = tic;
%test
Color=uint8(Color);
Color_=double(rgb2gray(Color));
RScolor = GetRelSmooth(Color_, rRelSmooth, rRelSmooth);
%over

% %test
figure(1);imshow(RScolor);

%depth图idea3
DepthUpdate=DepthUpdate/255;
r=1;
Sx= [diff(DepthUpdate,1,2), DepthUpdate(:,1,:) - DepthUpdate(:,end,:)];
Sy = [diff(DepthUpdate,1,1); DepthUpdate(1,:,:) - DepthUpdate(end,:,:)];
N = boxfilter(ones(m,n), r);
mean_Sx = boxfilter(Sx, r) ./ N;
mean_Sy = boxfilter(Sy, r) ./ N;
gradient_DepthUpdate=(sqrt(mean_Sx.^2+mean_Sy.^2)).*40;
 figure(12);imshow(gradient_DepthUpdate);
 
BW2 = ones(m, n);
BW2(gradient_DepthUpdate>0.07)=0;
figure(435);imshow(BW2);
BW2 = ones(m, n);
BW2(gradient_DepthUpdate>0.2)=0;
figure(436);imshow(BW2);
BW2 = ones(m, n);
BW2(gradient_DepthUpdate>0.15)=0;
figure(437);imshow(BW2);
%depth图idea3over
% patch1=BW2(841:1080,141:540);
% figure(420);imshow(patch1);
% patch2=BW2(141:380,840:1240);
% figure(421);imshow(patch2);

%color enhance测试
% Color1 = double(imread('E:\images\GFenhance_laundry.png'));
% Color1=uint8(Color1);
% Color1_=double(rgb2gray(Color1));
% RScolor1 = GetRelSmooth(Color1_, rRelSmooth, rRelSmooth);
% figure(10);imshow(RScolor1);
% Color1_=Color1_/255;
%     %color enhance idea3
% r=1;      
% Sx= [diff(Color1_,1,2), Color1_(:,1,:) - Color1_(:,end,:)];
% Sy = [diff(Color1_,1,1); Color1_(1,:,:) - Color1_(end,:,:)];
% N = boxfilter(ones(m,n), r);
% mean_Sx = boxfilter(Sx, r) ./ N;
% mean_Sy = boxfilter(Sy, r) ./ N;
% % figure(1111);imshow(mean_Sx.*40);
% % figure(1112);imshow(mean_Sy.*40);
% gradient_Color_=(sqrt(mean_Sx.^2+mean_Sy.^2)).*10;
% % figure(13);imshow(gradient_Color_);
% BW2 = ones(m, n);
% BW2(gradient_Color_>0.56)=0;
% figure(442);imshow(BW2);
% 
% BW = ones(m, n);
% BW((gradient_DepthUpdate<0.26)&(gradient_Color_>0.56)) = 0;
% figure(337);imshow(BW);
% %colorenhance over
% 
% % figure(2);imshow((Color_/255));
% % BW1 = ones(m, n);
% % BW1(RScolor<0.8) = 0;
% % figure(11);imshow(BW1);
% % %%测试over

RS = GetRelSmooth(DepthUpdate, rRelSmooth, rRelSmooth);
% 
%color图idea3
Color_=Color_/255;
r=1;
Sx= [diff(Color_,1,2), Color_(:,1,:) - Color_(:,end,:)];
Sy = [diff(Color_,1,1); Color_(1,:,:) - Color_(end,:,:)];
N = boxfilter(ones(m,n), r);
mean_Sx = boxfilter(Sx, r) ./ N;
mean_Sy = boxfilter(Sy, r) ./ N;
% % figure(1111);imshow(mean_Sx.*40);
% % figure(1112);imshow(mean_Sy.*40);
gradient_Color_=(sqrt(mean_Sx.^2+mean_Sy.^2)).*10;
% % figure(22);imshow(gradient_Color_);
% BW2 = ones(m, n);
% BW2(gradient_Color_>0.6)=0;
% figure(559);imshow(BW2);
% 
% % BW = ones(m, n);
% % BW((gradient_DepthUpdate<0.23)&(gradient_Color_>0.05)) = 0;
% % figure(339);imshow(BW);
%     %%color平滑，depth边缘
BW = ones(m, n);
BW((gradient_DepthUpdate<=0.15)&(gradient_Color_>=0.5)) = 0;
figure(33);imshow(BW);
BW = ones(m, n);
BW((gradient_DepthUpdate<=0.07)&(gradient_Color_>=0.5)) = 0;
figure(331);imshow(BW);
BW = ones(m, n);
BW((gradient_DepthUpdate>=0.38)&(gradient_Color_<=0.5)) = 0;
figure(332);imshow(BW);
% %color图idea3over
% 
% 
% BW = ones(m, n);
% % %测试专用
% % figure(3);imshow(RS);
% % BW((RS>0.96)&(RScolor<0.96)) = 0;
% % figure(3333);imshow(BW);
% % BW = ones(m, n);
% % BW((RS>0.96)&(RScolor<0.8)) = 0;
% % figure(333);imshow(BW);
% % %%color平滑，depth边缘
% % BW = ones(m, n);
% % BW((RS<0.96)&(RScolor>0.9)) = 0;
% % figure(34);imshow(BW);
% 
% %depth图idea3

BW2 = ones(m, n);
BW2(gradient_DepthUpdate>0.35)=0;
figure(446);imshow(BW2);
patch1=BW2(841:1080,141:540);
patch2=BW2(141:380,840:1240);
figure(420);imshow(patch1);
figure(421);imshow(patch2);
imwrite(patch1,'E:\ICASSIP2018\my paper\patchgradient\book\patch1_RS_35.png');
imwrite(patch2,'E:\ICASSIP2018\my paper\patchgradient\book\patch2_RS_35.png');
imwrite(BW2,'E:\ICASSIP2018\my paper\patchgradient\book\RS_35.png');

BW2 = ones(m, n);
BW2(gradient_DepthUpdate>0.19)=0;
figure(447);imshow(BW2);

BW2 = ones(m, n);
BW2(gradient_DepthUpdate>0.42)=0;
figure(448);imshow(BW2);
patch1=BW2(841:1080,141:540);
patch2=BW2(141:380,840:1240);
imwrite(patch1,'E:\ICASSIP2018\my paper\patchgradient\book\patch1_PG_38.png');
imwrite(patch2,'E:\ICASSIP2018\my paper\patchgradient\book\patch2_PG_38.png');
imwrite(BW2,'E:\ICASSIP2018\my paper\patchgradient\book\PG_38.png');

%depth图idea3over

BW2 = ones(m, n);
BW2(RS<0.88)=0;
figure(443);imshow(BW2);
patch1=BW2(841:1080,141:540);
patch2=BW2(141:380,840:1240);
imwrite(patch1,'E:\ICASSIP2018\my paper\patchgradient\book\patch1_RS_88.png');
imwrite(patch2,'E:\ICASSIP2018\my paper\patchgradient\book\patch2_RS_88.png');
imwrite(BW2,'E:\ICASSIP2018\my paper\patchgradient\book\RS_88.png');

BW2 = ones(m, n);
BW2(RS<0.92)=0;
figure(444);imshow(BW2);
% BW2 = ones(m, n);
% BW2(RS<0.98)=0;
% figure(445);imshow(BW2);
patch1=BW2(841:1080,141:540);
patch2=BW2(141:380,840:1240);
imwrite(patch1,'E:\ICASSIP2018\my paper\patchgradient\book\patch1_RS_92.png');
imwrite(patch2,'E:\ICASSIP2018\my paper\patchgradient\book\patch2_RS_92.png');
imwrite(BW2,'E:\ICASSIP2018\my paper\patchgradient\book\RS_92.png');

BW2 = ones(m, n);
BW2(RS<0.94)=0;
figure(440);imshow(BW2);
patch1=BW2(841:1080,141:540);
patch2=BW2(141:380,840:1240);
imwrite(patch1,'E:\ICASSIP2018\my paper\patchgradient\book\patch1_RS_94.png');
imwrite(patch2,'E:\ICASSIP2018\my paper\patchgradient\book\patch2_RS_94.png');
imwrite(BW2,'E:\ICASSIP2018\my paper\patchgradient\book\RS_94.png');

patch=BW2(140:540,840:1075);
figure(440);imshow(patch);

%%测试over
% BWsobel=edge(DepthUpdate,'sobel',2);
% figure(5);imshow(BWsobel);


BW(RS>=0.96) = BWAdp(1);
BW(RS<0.96) = BWAdp(2);
BW(RS<0.8) = BWAdp(3);
BW(RS<0.7) = BWAdp(4);
t12 = toc(t11);
fprintf('Computing relative smoothness costs %f s\n', t12);


%%%%%%%%%%% prepear LUT %%%%%%%%%%%%%
t11 = tic;
colorMax = 255; % both the color images and depth maps are normalized into [0, 255]
depthMax = 255;

%原来代码
% ColorRange = 0:3*(colorMax + 10)^2;
% ColorRange = ColorRange';

%自加：
ColorRange = 0:(colorMax + 10)^2;
ColorRange = ColorRange';
%over
DepthRange = 0:(depthMax + 10)^2;
DepthRange = DepthRange';

%original code
% ColorWeightLUT = exp(-ColorRange/(3*2*sigmaC^2));
%自加：
ColorWeightLUT = exp(-ColorRange/(2*sigmaC^2));
%over

%%测试专用
figure(1);
plot(ColorRange,ColorWeightLUT);
%%测试over

DepthWeightLUT = zeros(length(DepthRange), length(BWAdp));
for i = 1: length(BWAdp) 
    DepthWeightLUT(:, i) = exp(-DepthRange/(2*BWAdp(i)^2));
end
%%测试专用
figure(2);
plot(DepthRange,DepthWeightLUT(:,1));
figure(3);
plot(DepthRange,DepthWeightLUT(:,2));
figure(4);
plot(DepthRange,DepthWeightLUT(:,3));
figure(5);
plot(DepthRange,DepthWeightLUT(:,4));
%%测试over

t12 = toc(t11);
fprintf('Prepearing LUT costs %f s\n', t12);


%%%%% Color Weight %%%%%%%%
t11 = tic;
[ColorWeight] = mexGetColorWeight(Color(:,:,1), Color(:,:,2), Color(:,:,3), sigmaC, rSmooth, ColorWeightLUT);
t12=toc(t11);
fprintf('Computing color weight costs %f s\n', t12);


%%%%%%%%%%%

for i = 1: IterNumFixedPoint

    t11 = tic;
    [DepthWeightSmooth] = mexGetDepthWeight(DepthUpdate, rSmooth, BW, sigmaSSmooth, BWAdp, DepthWeightLUT);
    t12=toc(t11);
    fprintf('Computing depth weight costs %f s\n', t12);

    t11 = tic;
%原
%     [WeightSumSmooth, WeightedDepthUpd] = mexGetWeightedDepth(DepthUpdate, DepthWeightSmooth, ColorWeight, m, n, rSmooth);
%自加
     [WeightSumSmooth, WeightedDepthUpd] = mexGetWeightedDepth(DepthUpdate, DepthWeightSmooth, ColorWeight, m, n, rSmooth,RScolor,RS);
%over
     t12=toc(t11);
    fprintf('Computing weighted sum costs %f s\n', t12);


    t11 = tic;
    ResFixedPoint = ((1-alpha)*DepthNoise + 2*alpha*WeightedDepthUpd)./((1-alpha) + 2*alpha*WeightSumSmooth);
    t12=toc(t11);
    fprintf('Computing results using fixed point equation costs %f s\n', t12);

    DepthUpdate = ResFixedPoint;

    Diff = abs(Depth - DepthUpdate);
    MAE = sum(Diff(:))/(m*n);
    fprintf('Iteration %d, MAE is %f \n\n', i, MAE);


%     figure
    imshow(uint8(ResFixedPoint))

end

Res = ResFixedPoint;

