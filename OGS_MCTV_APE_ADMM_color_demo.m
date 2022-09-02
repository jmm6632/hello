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
 I =imread('peppers.tiff');
%  I =imread('4.1.05.tiff'); % house 


% OrgSigma= 30;

in=double(I);

% r通道
ur    =   in(:,:,1);
Nr     =   numel(ur);             
[m_r,n_r] =   size(ur);

%blur kernel definition
%  K     =   fspecial('average',15); %for denoising

K     =   fspecial('Gaussian',9,5); % for debluring
%  K   =   fspecial('motion',15,5);

blur_im_r  = imfilter(ur,K,'circular','conv');

BSNR_r= 40;%20*log10(norm(blur_im_r(:)-mean(blur_im_r(:)),'fro')/sqrt(Nr)/OrgSigma);
fprintf('BSNR of the observed image: %g dB.\n', BSNR_r);
  OrgSigma = BSNR2WGNsigma(blur_im_r, BSNR_r);

F_r        = blur_im_r + OrgSigma*randn(m_r,n_r);         %add noise


tao0=0.006; %for Gaussian/average blur this one is ok; for moving blur it should be larger
tao1=0.03;
% slightly tuning tao0 or tao1 may cause more appealing result 
if size(K)==1
tao  = -BSNR_r*tao1+1.09; %for denoising
else
tao  = -BSNR_r*tao0+1.09; %for deblurring
end
c    =  tao*m_r*n_r*OrgSigma.^2; % upper bound for the constraint


Param_r.OrigIm     = ur;      Param_r.MaxIter    = 1500; 
Param_r.SolRE      = 1e-5 ;Param_r.UpBound    = c;
Param_r.Beta       =2;       Param_r.Gamma      =1;
Param_r.Tao        = 1;       Param_r.BSNR       = BSNR_r;

%%***********************************************************
output_r = OGS_MCTV_APE_ADMM(F_r, K, Param_r); %% main program
%%***********************************************************
u_r     = output_r.Sol;           Reglambda  = output_r.Reglambda;
PSNR_r      = output_r.PSNR;     mse        = output_r.MSE;
IterTime_r   = output_r.IterTime; Fvalue     = output_r.Fvalue;

% g 通道
ug    =   in(:,:,2);
Ng     =   numel(ug);             
[m_g,n_g] =   size(ug);

%blur kernel definition
%  K     =   fspecial('average',15); %for denoising

% K     =   fspecial('Gaussian',9,5); % for debluring
%  K   =   fspecial('motion',15,5);

blur_im_g  = imfilter(ug,K,'circular','conv');

BSNR_g =40 ;%20*log10(norm(blur_im_g(:)-mean(blur_im_g(:)),'fro')/sqrt(Ng)/OrgSigma);

OrgSigma = BSNR2WGNsigma(blur_im_g, BSNR_g);

F_g        = blur_im_g + OrgSigma*randn(m_g,n_g);         %add noise


tao0=0.006; %for Gaussian/average blur this one is ok; for moving blur it should be larger
tao1=0.03;
% slightly tuning tao0 or tao1 may cause more appealing result 
if size(K)==1
tao  = -BSNR_g*tao1+1.09; %for denoising
else
tao  = -BSNR_g*tao0+1.09; %for deblurring
end
c    =  tao*m_r*n_r*OrgSigma.^2; % upper bound for the constraint


Param_g.OrigIm     = ug;       Param_g.MaxIter    = 1000; 
Param_g.SolRE      = 1e-5 ;    Param_g.UpBound    = c;
Param_g.Beta       =2;       Param_g.Gamma      =1;
Param_g.Tao        = 1;        Param_g.BSNR       = BSNR_g;

%%***********************************************************
output_g = OGS_MCTV_APE_ADMM(F_g, K, Param_g); %% main program
%%***********************************************************
u_g    = output_g.Sol;           Reglambda  = output_g.Reglambda;
PSNR_g       = output_g.PSNR;     mse        = output_g.MSE;
IterTime_g   = output_g.IterTime; Fvalue     = output_g.Fvalue;


% b通道
ub    =   in(:,:,3);
Nb     =   numel(ub);             
[m_b,n_b] =   size(ub);

%blur kernel definition


% K     =   fspecial('Gaussian',9,5); % for debluring
%  K   =   fspecial('motion',15,5);

blur_im_b  = imfilter(ub,K,'circular','conv');

BSNR_b =40;%20*log10(norm(blur_im_b(:)-mean(blur_im_b(:)),'fro')/sqrt(Nb)/OrgSigma);

OrgSigma = BSNR2WGNsigma(blur_im_b, BSNR_b);

F_b        = blur_im_b + OrgSigma*randn(m_b,n_b);         %add noise


tao0=0.006; %for Gaussian/average blur this one is ok; for moving blur it should be larger
tao1=0.03;
% slightly tuning tao0 or tao1 may cause more appealing result 
if size(K)==1
tao  = -BSNR_b*tao1+1.09; %for denoising
else
tao  = -BSNR_b*tao0+1.09; %for deblurring
end
c    =  tao*m_r*n_r*OrgSigma.^2; % upper bound for the constraint


Param_b.OrigIm     = ub;       Param_b.MaxIter    = 1500; 
Param_b.SolRE      = 1e-5 ;    Param_b.UpBound    = c;
Param_b.Beta       =2;       Param_b.Gamma      =1;
Param_b.Tao        = 1;        Param_b.BSNR       = BSNR_b;

%%***********************************************************
output_b = OGS_MCTV_APE_ADMM(F_b, K, Param_b); %% main program
%%***********************************************************
u_b    = output_b.Sol;           Reglambda  = output_b.Reglambda;
PSNR_b       = output_b.PSNR;     mse        = output_b.MSE;
IterTime_b  = output_b.IterTime; Fvalue     = output_b.Fvalue;



subplot(231);
imshow(I);title(sprintf('original '));


subplot(232);
% F=cat(3, F_r, F_g, F_b);
F_out=cat(3, F_r,F_g,F_b);
imshow(uint8(F_out),[]);
MSE_F=(MSE(F_r,ur)+MSE(F_g,ug)+MSE(F_b,ub))/3;
psnr_fun_F=(psnr_fun(F_r,ur)+psnr_fun(F_g,ug)+psnr_fun(F_b,ub))/3;
ssim_index_F=(ssim_index(F_r,ur)+ssim_index(F_g,ug)+ssim_index(F_b,ub))/3;

title(sprintf('noisy (MSE = %3.3f,PSNR = %3.3f dB,SSIM = %3.3f,C-score=%.3f)',MSE_F,psnr_fun_F,ssim_index_F,psnr_fun_F*ssim_index_F/MSE_F));
fprintf('noisy (MSE = %3.3f,PSNR = %3.3f dB,SSIM = %3.3f,C-score=%.3f)',MSE_F,psnr_fun_F,ssim_index_F,psnr_fun_F*ssim_index_F/MSE_F);

% imshow(uint8(F_r),[]);
% title(sprintf('noisy (MSE = %3.3f,PSNR = %3.3f dB,SSIM = %3.3f)',MSE(F_r,ur),psnr_fun(F_r,ur),ssim_index(F_r,ur)));
% fprintf('noisy (MSE = %3.3f,PSNR = %3.3f dB,SSIM = %3.3f)',MSE(F_r,ur),psnr_fun(F_r,ur),ssim_index(F_r,ur));

subplot(234);
 u=cat(3, u_r, u_g, u_b);
 imshow(uint8(u),[]);
MSE_u=(MSE(u_r,ur)+MSE(u_g,ug)+MSE(u_b,ub))/3;
psnr_fun_u=(psnr_fun(u_r,ur)+psnr_fun(u_g,ug)+psnr_fun(u_b,ub))/3;
ssim_index_u=(ssim_index(u_r,ur)+ssim_index(u_g,ug)+ssim_index(u_b,ub))/3;

title(sprintf('OGS_MCTV (MSE = %3.3f,PSNR = %3.3f dB,SSIM = %3.3f,C-score=%.3f,cputime= %.3f s)',MSE_u,psnr_fun_u,ssim_index_u,psnr_fun_u*ssim_index_u/MSE_u,max(IterTime_r)+max(IterTime_g)+max(IterTime_b)));
fprintf('OGS_MCTV (MSE = %3.3f,PSNR = %3.3f dB,SSIM = %3.3f,C-score=%.3f,cputime= %.3f s)',MSE_u,psnr_fun_u,ssim_index_u,psnr_fun_u*ssim_index_u/MSE_u,max(IterTime_r)+max(IterTime_g)+max(IterTime_b));
% imshow(uint8(u_r),[]);
% title(sprintf('OGS_MCTV (MSE = %3.3f,PSNR = %3.3f dB,SSIM = %3.3f,C-score= %.3f,cputime= %.3f s)',MSE(u_r,ur),psnr_fun(u_r,ur),ssim_index(u_r,ur),psnr_fun(u_r,ur)*ssim_index(u_r,ur)/MSE(u_r,ur),max(IterTime)));
% fprintf('OGS_MCTV (MSE = %3.3f,PSNR = %3.3f dB,SSIM = %3.3f,C-score= %.3f,cputime= %.3f s)',MSE(u_r,ur),psnr_fun(u_r,ur),ssim_index(u_r,ur),psnr_fun(u_r,ur)*ssim_index(u_r,ur)/MSE(u_r,ur),max(IterTime));

% title(sprintf('APE_ADMM (MSE = %3.3f,PSNR = %3.3f dB,SSIM = %3.3f,EBCM = %3.3f,C-score= %.3f,cputime= %.3f s)',MSE(u,double(I)),psnr_fun(u,double(I)),ssim_index(u,double(I)),EBCM(u),psnr_fun(u,double(I))*ssim_index(u,double(I))*EBCM(u)/MSE(u,double(I)),max(IterTime)));
% fprintf('APE_ADMM (MSE = %3.3f,PSNR = %3.3f dB,SSIM = %3.3f,EBCM = %3.3f,C-score= %.3f,cputime= %.3f s)',MSE(u,double(I)),psnr_fun(u,double(I)),ssim_index(u,double(I)),EBCM(u),psnr_fun(u,double(I))*ssim_index(u,double(I))*EBCM(u)/MSE(u,double(I)),max(IterTime));
% imwrite(uint8(u),'ogs_mctv_denoisy.png')

% 
% subplot(234);
%  plot(IterTime,PSNR,'-k.','LineWidth',2.5);
% 
% subplot(235);
% semilogy(IterTime,Reglambda,'--k*');

