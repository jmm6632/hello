function OutPut = OGS_MCTV_APE_ADMM(f, K, Param)
%  11/20/2021   by Ji Meimei 
%%
% ADMM method is applied to constrained OGS-MCTV-regularized image deconvolution problem.
% Suppose the mathematical model is given by
%
% $$f = K * u + n$$
% where f, u and n are the observed image with size $n\times m$, the original image 
% and the noise respectively. K is the blur matrix. 
%
% To restore f, We solve the constrained minimization problem
% 
% $$\min_u TV(u)$$
% subject to $||K*u-f||_2 \leq c$ .
%% Input:
%         Param.OrigIm:     original image                   
%         Param.MaxIter:    default=500;                
%         Param.Sol:        initial solution, default=f 
%         Param.Disp:       display the result or not. default=1;
%         Param.UpBound:    upper bound for the norm of the residual $c$ 
%         Param.Reglambda:  regularization parameter  $\lambda$
%         Param.SolRE:      stop criterion for the relative difference.   default=1e-4.
%         Param.Beta        ADMM parameter.   default=1.
%         Param.Gamma       ADMM parameter.   default=1.618
%%  OutPut:
%         OutPut.Sol:       restored image u;
%         OutPut.PSNR:      the evolution of PSNR against iteration number
%         OutPut.IterTime:  the evolution of CPU time against iteration
%         number

%%=========================================================================
%  Copyright(c), Jun. 2013, Dr.HE Chuan(hechuan8512@163.com)
%%========================================================================= 
u       = f;               MaxIter = 500;       
SolRE   = 1e-5;            beta    = 1;          gamma = 1.618;
gcv     = 0;               flag    = 1;          tao   = 1;
[m,n]   = size(f);
kesai1  = zeros(m,n);      kesai2  = kesai1;     mu    = kesai1;
deta    = mu  ;
z = zeros(m,n);
if nargin == 3
    if isfield(Param,'OrigIm'),     xtrue     = Param.OrigIm;     end
    if isfield(Param,'MaxIter'),    MaxIter   = Param.MaxIter;    end
    if isfield(Param,'Beta'),       beta      = Param.Beta;       end                
    if isfield(Param,'Gamma'),      gamma     = Param.Gamma;      end
    if isfield(Param,'Tao'),        tao       = Param.Tao;        end
    if isfield(Param,'BSNR'),       BSNR      = Param.BSNR;       end
    if isfield(Param,'Disp'),       flag      = Param.Disp;       end
    if isfield(Param,'SolRE'),      SolRE     = Param.SolRE;      end
    if isfield(Param,'Reglambda'),  lambda    = Param.Reglambda;  end
    if isfield(Param,'UpBound'),    c         = Param.UpBound;    
                                    RegMethod = 'Discrepancy';    end
    if isfield(Param,'GCV'),        gcv       = Param.GCV;        
                                    RegMethod = 'GCV';            end
    if isfield(Param,'Sol'),        u         = Param.Sol;        end
end

beta1     = 1* beta;% 10^(BSNR/10-1)*beta;
beta2     =  1* beta ; %10^(BSNR/10-1)*beta;
beta3     = 1* beta;
 alpha     =0.06;
grpSz     =3;
Nit_inner =5;

if ~isvar('c') && ~isvar('lambda') && ~gcv 
     sigma     = ImageStdDev(f);%noise variance estimation
     c         = tao*sigma.^2 * m * n;
     RegMethod = 'Discrepancy';
end

if isvar('lambda')
     RegMethod = 'RegFixed';
end


Param.Sol  = u;
OutPut     = Param;

X = zeros(m,n);V1 = zeros(m,n); V2=V1;
C = getC;
[D,Dt] = defDDt;

regpar = zeros(MaxIter,1);   PSNR=regpar;  iter_time = regpar; mse=regpar;
Fvalue = zeros(MaxIter,1);
cont   = 1; k = 0;  
%finite diff
[D1U,D2U] = D(u);

while cont
    tic;
    k = k + 1;
     
     % ==================
     % Y-subprolem Shrinkage Step
     % ==================
     Y1 = gstvdm(beta2/(beta2-alpha)*D1U + (kesai1 - alpha*V1)/(beta2-alpha) , grpSz , 1/(beta2-alpha), Nit_inner);
     Y2 = gstvdm(beta2/(beta2-alpha)*D2U + (kesai2 - alpha*V2)/(beta2-alpha) , grpSz , 1/(beta2-alpha), Nit_inner);
     % ==================
     % V-subprolem Shrinkage Step
     V1 = gstvdm(Y1 , grpSz , 1/alpha, Nit_inner);
     V2 = gstvdm(Y2 , grpSz , 1/alpha, Nit_inner);
    % ==================
    %     z-subprolem
    % ==================
        z = min(255,max(u + deta/beta3,0));
    % ==================
    %     U-subprolem
    % ==================
    KtF  = (beta1/beta2)*(conj(C.eigsK) .*fft2(X-mu/beta1));
    unew = KtF + fft2(Dt(Y1-kesai1/beta2,Y2-kesai2/beta2) + (beta1/beta2)*(z - deta));
    Unew = unew./(C.eigsDtD + (beta1/beta2)*C.eigsKtK + beta3/beta2);
    unew = real(ifft2(Unew));
    
    % ==================
    %   X-subprolem 
    % ==================
    switch RegMethod
        case 'Discrepancy'
            if  size(K)==1 %for denoising
                Kunew   =  real(ifft2(Unew));  
            else  % for deblurring
                Kunew   =  real(ifft2(C.eigsK .*Unew));
            end
            
            A       =  Kunew + mu/beta1;
            FA_norm =  norm(f-A,'fro');
            
        if  FA_norm^2<=c
            lambda  =  0;
            X       =  A;
        else
            lambda  =  beta1* FA_norm/sqrt(c)-beta1; 
            X       =  (lambda*f + beta1*A)/(lambda + beta1);
            
            
        end
            
            
        case 'RegFixed'
            Kunew   =  real(ifft2(C.eigsK .*Unew));
            A       =  Kunew + mu/beta1;
            X       = (lambda*f + beta1*A)/(lambda + beta1);
                 
    end
            
    %% update kesai an mu
    [D1U,D2U] = D(unew);
     mu       = mu     - gamma*beta1*(X  - Kunew);
     kesai1   = kesai1 - gamma*beta2*(Y1 - D1U);
     kesai2   = kesai2 - gamma*beta2*(Y2 - D2U);  
     
     deta    =  deta -gamma* beta3* (z - unew) ;
    %% relative error and checking stopping rule 
    re   = norm(unew-u,'fro')/norm(u,'fro');
%      if re<1e-4
%        beta1=1.2*beta1;
%        beta2=1.2*beta2;% slightly tuning may accelerate the convergence 
%     end
    u    = unew;
    fvalue = fval;
    
    cont = (k<MaxIter)&&(re>SolRE);    
    if isvar('xtrue'),     
        PSNR(k) = psnr(u,xtrue);
        if mod(k,10)==0 && flag
            fprintf('%3d-th  PSNR: %2.2f,   regpara: %1.2e  re: %1.2e\n', k,PSNR(k), lambda, re);
        end
    end
    tElapsed     = toc;
    if  k==1       
         iter_time(k)  =tElapsed; 
    else iter_time(k) = iter_time(k-1)+tElapsed;
    end
    Fvalue(k)    = fvalue;
    regpar(k)    = lambda;
    mse(k)       = norm(unew-xtrue,'fro')^2/m/n;
end
OutPut.Sol       = u;
OutPut.PSNR      = PSNR(1:k);
OutPut.IterTime  = iter_time(1:k);
OutPut.Reglambda = regpar(1:k);
OutPut.Fvalue    = Fvalue(1:k);
OutPut.MSE       = mse(1:k);
            
%% function used above 
function C = getC
        sizeF = size(f);
        C.eigsK = psf2otf(K,sizeF);
        C.eigsDtD1 =  abs(psf2otf([1,-1],sizeF)).^2;
        C.eigsDtD2 =  abs(psf2otf([1;-1],sizeF)).^2;
        C.eigsDtD  =  C.eigsDtD1 + C.eigsDtD2;
        C.eigsKtK = abs(C.eigsK).^2;
end



function Fvalue = fval
        Fvalue = sum(sum(sqrt(D1U.^2 + D2U.^2)));
        KU_F = real(ifft2(C.eigsK .* fft2(u))) - f;
        Fvalue = Fvalue + lambda/2 * norm(KU_F,'fro')^2;
end


function [D,Dt] = defDDt
        % defines finite difference operator D
        % and its transpose operator
        % referring to FTVD code
        
        D = @(U) ForwardD(U);
        Dt = @(X,Y) Dive(X,Y);
end
        
function [Dux,Duy] = ForwardD(U)
         % Forward finite difference operator
         Dux = [diff(U,1,2), U(:,1) - U(:,end)];
         Duy = [diff(U,1,1); U(1,:) - U(end,:)];
end
        
function DtXY = Dive(X,Y)  %Dt=-div
        % Transpose of the forward finite difference operator
        DtXY = [X(:,end) - X(:, 1), -diff(X,1,2)];
        DtXY = DtXY + [Y(end,:) - Y(1, :); -diff(Y,1,1)];
end

function sigma = ImageStdDev(y)
h = [0.03522629188571 0.08544127388203 -0.13501102001025 -0.45987750211849 0.80689150931109 -0.33267055295008];       
h = h(end:-1:1); %% flip is used only to give exactly same result as previous version of the code (which was using filter2)
        
z = conv2(y,h,'same');
z=conv2(z,h','same');
sigma = median(abs(z(:)))/.6745;
end

function tf = isvar(name)
%function tf = isvar(name)
% determine if "name" is a variable in the caller's workspace

if nargin < 1
	help isvar
	error arg
end

tf = true(1);
evalin('caller', [name ';'], 'tf=logical(0);')
end

            
            
            
            
 
    
    

            
            
 

end