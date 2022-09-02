%demo of APE_ADMM algorithm (TV regualrized ADMM with adaptive parameter estimation)

clc; clear all; close all;
path(path,genpath(pwd));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating the observed image
 
%  I     =   imread('Einstein256.bmp');
%   I =imresize(imread('peppers.tiff'),[256,256]);
%  I =imread('Cameraman256.png');
%  I =imread('parrots.tif');
%  I     =   imread('Lena.bmp');
%  I     =   imread('butterfly2.bmp');
%  I     =   imread('lady_liberty.png');
% I     =   imread('shroom.png');
%   I =imread('bridge.png');
   I     =   imread('snow_leaves.png');
%    I =imread('child_swimming.png'); 
%  I     =   imread('goldhill.bmp');
%  I     =   imread('boat512.tiff');
%  I     =   imread('man.tiff');
%  OrgSigma= 30;

if size(I,3) > 1
    I = rgb2gray(I);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
end
u0    =   double(I);
N     =   numel(u0);             
[m,n] =   size(u0);

%blur kernel definition
%  K     =   fspecial('average',15); %for denoising

K     =   fspecial('Gaussian',9,5); % for debluring
%  K   =   fspecial('motion',15,5);

blur_im  = imfilter(u0,K,'circular','conv');

BSNR=40 ;%20*log10(norm(blur_im(:)-mean(blur_im(:)),'fro')/sqrt(N)/OrgSigma);
fprintf('BSNR of the observed image: %g dB.\n', BSNR);
OrgSigma = BSNR2WGNsigma(blur_im, BSNR);

F        = blur_im + OrgSigma*randn(m,n);         %add noise
PSNR_F    =psnr(F,u0);
fprintf('PSNR of the observed image: %g dB.\n\n', PSNR_F);


tao0=0.006; %for Gaussian/average blur this one is ok; for moving blur it should be larger
tao1=0.05;
% slightly tuning tao0 or tao1 may cause more appealing result 
if size(K)==1
tao  = -BSNR*tao1+1.09; %for denoising
else
tao  = -BSNR*tao0+1.09; %for deblurring
end
c    =  tao*m*n*OrgSigma.^2; % upper bound for the constraint


Param.OrigIm     = u0;      Param.MaxIter    =1000; 
Param.SolRE      = 1e-5 ;Param.UpBound    = c;
Param.Beta       =2;       Param.Gamma      =1;
Param.Tao        = 1;       Param.BSNR       = BSNR;

%%***********************************************************
output = OGS_MCTV_APE_ADMM(F, K, Param); %% main program
%%***********************************************************
u     = output.Sol;           Reglambda  = output.Reglambda;
PSNR       = output.PSNR;     mse        = output.MSE;
IterTime   = output.IterTime; Fvalue     = output.Fvalue;
% fprintf('Proposed APE_ADMM: Elapsed time is %g seconds.\n', IterTime(end));
% fprintf('Proposed APE_ADMM: Total iterations is %g.\n', length(IterTime));
% fprintf('Proposed APE_ADMM: Final regularization parameter is %g.\n', Reglambda(end));
% fprintf('Proposed APE_ADMM: Final PSNR is %g dB.\n\n\n', PSNR(end));

subplot(231);
imshow(I);title(sprintf('original '));

subplot(232);
imshow(uint8(F),[]);
title(sprintf('noisy (MSE = %3.3f,PSNR = %3.3f dB,SSIM = %3.3f)',MSE(F,double(I)),psnr_fun(F,double(I)),ssim_index(F,double(I))));
fprintf('noisy (MSE = %3.3f,PSNR = %3.3f dB,SSIM = %3.3f)',MSE(F,double(I)),psnr_fun(F,double(I)),ssim_index(F,double(I)));

subplot(234);
imshow(uint8(u),[]);
title(sprintf('OGS_MCTV (MSE = %3.3f,PSNR = %3.3f dB,SSIM = %3.3f,C-score= %.3f,cputime= %.3f s)',MSE(u,double(I)),psnr_fun(u,double(I)),ssim_index(u,double(I)),psnr_fun(u,double(I))*ssim_index(u,double(I))/MSE(u,double(I)),max(IterTime)));
fprintf('OGS_MCTV (MSE = %3.3f,PSNR = %3.3f dB,SSIM = %3.3f,C-score= %.3f,cputime= %.3f s)',MSE(u,double(I)),psnr_fun(u,double(I)),ssim_index(u,double(I)),psnr_fun(u,double(I))*ssim_index(u,double(I))/MSE(u,double(I)),max(IterTime));

% title(sprintf('APE_ADMM (MSE = %3.3f,PSNR = %3.3f dB,SSIM = %3.3f,EBCM = %3.3f,C-score= %.3f,cputime= %.3f s)',MSE(u,double(I)),psnr_fun(u,double(I)),ssim_index(u,double(I)),EBCM(u),psnr_fun(u,double(I))*ssim_index(u,double(I))*EBCM(u)/MSE(u,double(I)),max(IterTime)));
% fprintf('APE_ADMM (MSE = %3.3f,PSNR = %3.3f dB,SSIM = %3.3f,EBCM = %3.3f,C-score= %.3f,cputime= %.3f s)',MSE(u,double(I)),psnr_fun(u,double(I)),ssim_index(u,double(I)),EBCM(u),psnr_fun(u,double(I))*ssim_index(u,double(I))*EBCM(u)/MSE(u,double(I)),max(IterTime));
% imwrite(uint8(u),'ogs_mctv_denoisy.png')

% 
% subplot(234);
%  plot(IterTime,PSNR,'-k.','LineWidth',2.5);
% 
% subplot(235);
% semilogy(IterTime,Reglambda,'--k*');

