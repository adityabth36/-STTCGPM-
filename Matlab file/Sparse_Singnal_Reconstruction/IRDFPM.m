
%   F(x)=0. 
%
function [x1,mses,Tcpu,NF,NormF] = IRDFPM(A,b,lambda,true_x,method,model,para) 
 
format long

% if A is a matrix, we find out dimensions of y and x,
% and create function handles for multiplication by A and A',
% so that the code below doesn't have to distinguish between
% the handle/not-handle cases
if ~isa(A, 'function_handle')
  AT = @(x) (x'*A)'; %A'*x;
  A = @(x) A*x;
end
% from this point down, A and AT are always function handles.

% Precompute A'*b since it'll be used a lot
Atb = AT(b);    
x0 = Atb;
n = length(x0);

% initial point
xu0 =  x0.*(x0 >= 0);
xv0 = -x0.*(x0 <  0);
% from these two relations, x0 = u0-v0;

%% parameter
% mu =10;
% lambda_k=15;

%% the stopping criterion
epsilon = 1e-6;

%% the line search parameters and relaxation factor
k_max = para.Itr_max;   % the maximum number of iterations
gamma = para.gamma;     % the initial guess
sigma = para.sigma;     % the coefficient of line search 
tau = para.tau;         % the compression ratio
alpha = para.alpha;     % the coefficient of the inertial step
rho = para.rho;         % the relaxation factor 

fprintf('%s & LSmodel=%d & gamma=%.4f & sigma=%.4f & tau=%.4f & alpha=%.4f & rho=%.4f\n', ... 
    method,model,gamma,sigma,tau,alpha,rho);

% start the clock
t0 = cputime;

% % define function handle
% Fu = @(x,xu) min(xu,AT(A(x))-Atb+lambda);
% Fv = @(x,xv) min(xv,-AT(A(x))+Atb+lambda);

%% compute the search direction
Ax0 = A(x0);
tempx = AT(Ax0)-Atb;
Fxu0 = min(xu0,tempx+lambda);
Fxv0 = min(xv0,-tempx+lambda);
NF = 1;  
NormFxk2 = Fxu0'*Fxu0+Fxv0'*Fxv0;    
NormFxk = sqrt(NormFxk2);
xu0_old = xu0;
xv0_old = xv0;
L1 = 0;
     
for k=1:k_max
    
    if k==1 && NormFxk<=epsilon
        L1 = 1;
        NormF = NormFxk; % the final norm of equations
        mses(k) = 1/n*norm(x0-true_x)^2;
        Tcpu(k) = cputime-t0;
        break; 
    end

   
    %% compute the inertial step %%
    yuk = xu0+alpha*(xu0-xu0_old);
    yvk = xv0+alpha*(xv0-xv0_old);
    yk = yuk-yvk;
    Ayk = A(yk);
    tempy = AT(Ayk)-Atb;
    Fyuk = min(yuk,tempy+lambda);
    Fyvk = min(yvk,-tempy+lambda);
    NF = NF+1;
    NormFyk2 = Fyuk'*Fyuk+Fyvk'*Fyvk;
    NormFyk = sqrt(NormFyk2);
    if NormFyk<=epsilon
        L1 = 1;
        NormF = NormFyk;   % the final norm of equations
        mses(k) = 1/n*norm(yk-true_x)^2;
        Tcpu(k) = cputime-t0;
        break; 
    end
    
    %% compute the initial direction %%
    if k==1
        dku = -Fyuk;
        dkv = -Fyvk;
    else
        % update the search direction
        switch method
             case'FITTCGPM-PRP'
                 wuk = Fyuk-Fyuk_old;
                 wvk = Fyvk-Fyvk_old;
                 NormFyk2=Fyuk'*Fyuk+Fyvk'*Fyvk;
                 NormFyk=sqrt(NormFyk2);
                 NormFyk2_old = Fyuk_old'*Fyuk_old+Fyvk_old'*Fyvk_old;
                 NormFyk_old = sqrt(NormFyk2_old);
                 Normdk2=dku'*dku+dkv'*dkv;
                 Normdk=sqrt(Normdk2);
                 Fyktwk= Fyuk'*wuk+Fyvk'*wvk;
                 betak=Fyktwk/NormFyk2_old;
                 FyktFyk_old=Fyuk'*Fyuk_old+Fyvk'*Fyvk_old;
                 thetak=-0.055*FyktFyk_old/NormFyk2_old;%-0.34545
                 muk=0.001*NormFyk/max(NormFyk_old,(abs(betak))*Normdk);
                 dku=-Fyuk+muk*betak*dku+thetak*Fyuk_old;
                 dkv=-Fyvk+muk*betak*dkv+thetak*Fyvk_old;      
             case'FITTCGPM-DY'
                 wuk = Fyuk-Fyuk_old;
                 wvk = Fyvk-Fyvk_old;
                 NormFyk2=Fyuk'*Fyuk+Fyvk'*Fyvk;
                 NormFyk=sqrt(NormFyk2);
                 NormFyk2_old = Fyuk_old'*Fyuk_old+Fyvk_old'*Fyvk_old;
                 NormFyk_old = sqrt(NormFyk2_old);
                 Normdk2=dku'*dku+dkv'*dkv;
                 Normdk=sqrt(Normdk2);
                 dktwk=dku'*wuk+dkv'*wvk;
                 betak=NormFyk2_old/dktwk;
                 FyktFyk_old=Fyuk'*Fyuk_old+Fyvk'*Fyvk_old;
                 thetak=-0.055*FyktFyk_old/NormFyk2_old;%-0.34545
                 muk=0.001*NormFyk/max(NormFyk_old,abs(betak)*Normdk);
                 dku=-Fyuk+muk*betak*dku+thetak*Fyuk_old;
                 dkv=-Fyvk+muk*betak*dkv+thetak*Fyvk_old;     
            case'IMSMNE'
                t1=0.3; t2=10^-3; D=10; r0=3; tauk=6.2 ; 
                wuk = Fyuk-Fyuk_old;
                wvk = Fyvk-Fyvk_old;
                suk= yuk -yuk_old;
                svk=  yvk -yvk_old;
                normFyk_old2 = Fyuk_old(:)'*Fyuk_old(:)+Fyvk_old(:)'*Fyvk_old(:);
                normFyk_old = sqrt(normFyk_old2);
                normsk = suk(:)'*suk(:) +svk(:)'*svk(:);
                sktw0 = suk(:)'*wuk + svk(:)'*wvk(:) ;
                huk = wuk + D*normFyk_old^r0 *suk + max(0,((-sktw0)/(normsk)))*suk; 
                hvk = wvk + D*normFyk_old^r0 *svk + max(0,((-sktw0)/(normsk)))*svk; 
                normhk2 =huk(:)'*huk(:) + hvk(:)'*hvk(:) ; 
                normhk =sqrt(normhk2) ;
                Normdk2=dku(:)'*dku(:)+dkv(:)'*dkv(:);
                Normdk=sqrt(Normdk2);
                dkthk = dku(:)'*huk(:) +dkv(:)'*hvk(:);
                Fykthk = Fyuk(:)'*huk(:) +Fyvk(:)'*hvk(:);
                Fyktdk = Fyuk(:)'*dku(:) +Fyvk(:)'*dkv(:);
                k_a = max(tauk*normhk*Normdk ,dkthk);
                betak = Fykthk/k_a - ((normhk2 * Fyktdk)/(k_a^2));
                skthk = suk(:)'*huk(:) +svk(:)'*hvk(:) ;
                theta = skthk/normsk ;
                thetak= min(max(theta, t1), t2);
                dku = -thetak*Fyuk + betak* dku ;
                dkv = -thetak*Fyvk + betak* dkv ;
             case'STTCGPM' 
                 mu=10; lambda_k=15;
                 wuk = Fyuk-Fyuk_old;
                 wvk = Fyvk-Fyvk_old;
                 NormFyk2=Fyuk(:)'*Fyuk(:)+Fyvk(:)'*Fyvk(:);
                 % NormFyk=sqrt(NormFyk2);
                 Normdk2=dku(:)'*dku(:)+dkv(:)'*dkv(:);
                 Normdk=sqrt(Normdk2);
                 dktwk = dku(:)'*wuk(:)+dku(:)'*wvk(:);
                 Fyktdk = Fyuk(:)'*dku(:)+Fyvk(:)'*dkv(:);
                 Uk = mu *(Normdk2 + NormFyk2+ abs(dktwk));
                 betak = NormFyk2/Uk - ((NormFyk2 * (Fyktdk))/ Uk^2);
                 theta_k= lambda_k+betak*((Fyktdk)/(NormFyk2));
                 puk = Fyuk ;
                 pvk =Fyvk ;
                 Normpk2 = puk(:)'*puk(:) + pvk(:)'*pvk(:);
                 Normpk = sqrt(Normpk2);
                 dku=-theta_k*Fyuk+betak*dku+((Fyktdk)/(Normdk*Normpk))*puk;
                 dkv=-theta_k*Fyvk+betak*dkv+((Fyktdk)/(Normdk*Normpk))*pvk; 
        end
    end
    Normdk2 = dku'*dku+dkv'*dkv;
    Normdk = sqrt(Normdk2);
    
    
    %%% Start Armijo-type line search  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model=1 means -F(zk)'*dk    sigma*tk*norm(dk)^2
    % model=2 means -F(zk)'*dk    sigma*tk*norm(F(zk))*norm(dk)^2

    if model==1
        t = gamma;
        zuk_new = yuk+t*dku;
        zvk_new = yvk+t*dkv;
        zk_new = zuk_new-zvk_new;
        Azk_new = A(zk_new);
        tempz = AT(Azk_new)-Atb;
        Fzuk_new = min(zuk_new,tempz+lambda);
        Fzvk_new = min(zvk_new,-tempz+lambda);
        NF = NF+1;
        Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
        % check the Armijo-type line search condition
        while -Fzk_newtdk < sigma*t*Normdk2 && t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            zuk_new = yuk+t*dku;
            zvk_new = yvk+t*dkv;
            zk_new = zuk_new-zvk_new;
            Azk_new = A(zk_new);
            tempz = AT(Azk_new)-Atb;
            Fzuk_new = min(zuk_new,tempz+lambda);
            Fzvk_new = min(zvk_new,-tempz+lambda);
            NF = NF+1;
            Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
        end %%% End Armijo-type line search %%%
        NormFzk_new2 = Fzuk_new'*Fzuk_new+Fzvk_new'*Fzvk_new;
        NormFzk_new = sqrt(NormFzk_new2);
    else 
        t = gamma;
        zuk_new = yuk+t*dku;
        zvk_new = yvk+t*dkv;
        zk_new = zuk_new-zvk_new;
        Azk_new = A(zk_new);
        tempz = AT(Azk_new)-Atb;
        Fzuk_new = min(zuk_new,tempz+lambda);
        Fzvk_new = min(zvk_new,-tempz+lambda);
        NF = NF+1;
        Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
        NormFzk_new2 = Fzuk_new'*Fzuk_new+Fzvk_new'*Fzvk_new;
        NormFzk_new = sqrt(NormFzk_new2);
        % check the Armijo-type line search condition
        while -Fzk_newtdk < sigma*t*NormFzk_new*Normdk2 && t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            zuk_new = yuk+t*dku;
            zvk_new = yvk+t*dkv;
            zk_new = zuk_new-zvk_new;
            Azk_new = A(zk_new);
            tempz = AT(Azk_new)-Atb;
            Fzuk_new = min(zuk_new,tempz+lambda);
            Fzvk_new = min(zvk_new,-tempz+lambda);
            NF = NF+1;
            Fzk_newtdk = Fzuk_new'*dku+Fzvk_new'*dkv;
            NormFzk_new2 = Fzuk_new'*Fzuk_new+Fzvk_new'*Fzvk_new;
            NormFzk_new = sqrt(NormFzk_new2);
        end %%% End Armijo-type line search %%
    end
    zuk = zuk_new;
    zvk = zvk_new;
    zk = zk_new;
    Fzuk = Fzuk_new;
    Fzvk = Fzvk_new;
    NormFzk2 = NormFzk_new2;
    NormFzk = NormFzk_new;
    if NormFzk<=epsilon
        L1 = 1;
        NormF = NormFzk; % the final norm of equations
        mses(k) = 1/n*norm(zk-true_x)^2;
        Tcpu(k) = cputime-t0;
        break;
    end
    Fzktykzk = Fzuk'*(yuk-zuk)+Fzvk'*(yvk-zvk);
    xik = Fzktykzk/NormFzk2;
    % compute the next iteration 
    xu1 = yuk-rho*xik*Fzuk;
    xv1 = yvk-rho*xik*Fzvk;
%     xuv1min = min(xu1,xv1);
%     xu1 = xu1-xuv1min;
%     xv1 = xv1-xuv1min;
    x1 = xu1-xv1;
    mses(k) = 1/n*norm(x1-true_x)^2;
    Ax1 = A(x1);
    tempx = AT(Ax1)-Atb;
    Fxu = min(xu1,tempx+lambda);
    Fxv = min(xv1,-tempx+lambda);
    NF = NF+1;
    NormFx2 = Fxu'*Fxu+Fxv'*Fxv;
    NormFx = sqrt(NormFx2);
    if NormFx<=epsilon
        L1 = 1;
        NormF = NormFx;
        mses(k) = 1/n*norm(x1-true_x)^2;
        Tcpu(k) = cputime-t0;
        break;
    end
    
    % update the iteration
    xu0_old = xu0;
    xv0_old = xv0;
    xu0 = xu1;
    xv0 = xv1;
    yuk_old = yuk;
    yvk_old = yvk;
    Fyuk_old = Fyuk;
    Fyvk_old = Fyvk;
    NormFyk2_old = NormFyk2;
    Tcpu(k) = cputime-t0;
end
if L1~=1
    NF = NaN;
    NormF = NaN;
end
