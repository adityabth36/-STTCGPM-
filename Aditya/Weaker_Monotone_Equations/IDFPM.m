% Matlab Model by Jianghua Yin (Jan.,2022, Nanning)
% Copyright (C) 2022 Jian Group
% All Rights Reserved
%
%% the inertial derivative-free projection method (IDFPM) for solving 
%% constrained nonlinear pseudo-monotone equations of the form
%   F(x)=0, x\in C, 
% where C is a nonempty closed convex set.
%
function [Tcpu,NF,Itr,NormF] = IDFPM(NO,method,Switch,model,x0,para) 
 
format long

% start the clock
tic;

%% the number of itrations
% Itr=0;   

%% Initial information
[nprob,~]=init(NO);

%% the stopping criterion
epsilon=1e-6;
epsilon1=1e-7;

%% the line search parameters and relaxation factor
k_max = para.Itr_max;   % the maximum number of iterations
gamma = para.gamma;     % the initial guess
sigma = para.sigma;     % the coefficient of line search 
tau = para.tau;         % the compression ratio
 alpha = para.alpha;     % the coefficient of the inertial step
rho = para.rho;         % the relaxation factor 
 alpha_try = alpha; %0.1;
mu =10;
lambda_k=15;

 fprintf('%s & %s & LSmodel=%d & gamma=%.4f & sigma=%.4f & tau=%.4f & Switch=%.4f & rho=%.4f\n', ... 
     nprob,method,model,gamma,sigma,tau,Switch,rho);

%% compute the search direction
Fx0 = feval(nprob,x0,1);   % evaluate the function value specified by nprob at x0
NF = 1;  
NormFx0 = norm(Fx0);                   
x0_old = x0;
L1 = 0;
     
for k=1:k_max
    
    if k==1 && NormFx0<=epsilon
        L1 = 1;
        NormF = NormFx0; % the final norm of equations
        break; 
    end
    if k==1
        alpha = alpha_try;
    else
        if Switch==1
            alpha = min(alpha_try,1/(norm(x0-x0_old)*k^2));
        elseif Switch==2
            alpha = min(alpha_try,1/(k*norm(x0-x0_old))^2);
        else
            alpha = alpha_try;
        end
    end
    %% compute the inertial step %%
    v0 = x0+alpha*(x0-x0_old);
    Fv0 = feval(nprob,v0,1);
    NF = NF+1;
    NormFv0 = norm(Fv0);
    if NormFv0<=epsilon
        L1 = 1;
        NormF = NormFv0;   % the final norm of equations
        break; 
    end
    
    %% compute the initial direction %%
    if k==1
        dk = -Fv0;
    else
        % update the search direction
        switch method   
            case 'STTCGPM'
               w0 = Fv0-Fy0_old;
               U_k = mu *(norm(d_k_prev)^2 + NormFv0^2+ abs(d_k_prev'*w0));
               beta_k = NormFv0^2/U_k - ((NormFv0^2 * (Fv0'*d_k_prev))/ U_k^2);
               theta_k= lambda_k+beta_k*((Fv0'*d_k_prev)/(NormFv0^2));
               p_k = Fv0 ;
               dk=-theta_k.*Fv0 + beta_k .*d_k_prev + ((Fv0'*d_k_prev)/(norm(d_k_prev)*norm(p_k)))*p_k;
           case 'FITTCGPM-PRP'%第二篇
                 w0 = Fy0-Fy0_old;
                 betak=Fy0'*w0/(norm(Fy0_old)^2);
                 nuk=norm(Fy0)/max(norm(Fy0_old),(abs(betak))*norm(dk));
                 thetak=-0.055*Fy0'*Fy0_old/norm(Fy0_old)^2;  % 原始的分母是 norm(Fk0)^2  
                 dk=-Fy0+0.001*nuk*betak*dk+thetak*Fy0_old;
           case'FITTCGPM-DY'%第二篇
                 w0 = Fy0-Fy0_old;
                 betak=norm(Fy0)^2/(dk'*w0);
                 nuk=norm(Fy0)/max(norm(Fy0_old),(abs(betak))*norm(dk));
                 thetak=-0.055*Fy0'*Fy0_old/norm(Fy0_old)^2;  % 原始的分母是 norm(Fk0)^2  
                 dk=-Fy0+0.001*nuk*betak*dk+thetak*Fy0_old; 
%            case 'FITTCGPM-FR'%第二篇
%                 betak=norm(Fy0)^2/norm(Fy0_old)^2;
%                 nuk=norm(Fy0)/max(norm(Fy0_old),(abs(betak))*norm(dk));
%                 thetak=-0.055*Fy0'*Fy0_old/norm(Fy0_old)^2;  % 原始的分母是 norm(Fk0)^2  
%                 dk=-Fy0+0.001*nuk*betak*dk+thetak*Fy0_old;    
          case 'GITDFPA'
                s_k =v0-v0_old ;
                w0 = Fv0-Fy0_old;
                p_k = Fv0 ;
               U_k = 2.5 *(norm(d_k_prev)^2 + NormFv0^2)+ max([norm(d_k_prev)^2,norm(Fy0_old)^2, d_k_prev'*w0] );
             
               W_k = max(0.2*norm(d_k_prev)*norm(p_k), norm(Fy0_old)^2);
               eta_k = min(0.1 , max(0, (p_k'*(w0 -s_k))/(norm(p_k)^2)));
               beta_k = NormFv0^2/U_k - ((NormFv0^2 * (Fv0'*d_k_prev))/ U_k^2);
              
               dk=-Fv0 + beta_k .*d_k_prev +eta_k* ((Fv0'*d_k_prev)/(W_k))*p_k;
            otherwise
                disp('Input error! Please check the input method');
        end
    end
    Normdk = norm(dk);
    if Normdk<epsilon1
        L1 = 1;
        NormF = NormFy0;
        break;
    end
    Normdk2 = Normdk^2;
    Fy0_old = Fv0;
    v0_old =v0;
    %%% Start Armijo-type line search  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model=1 means -F(zk)'*dk ≥ sigma*tk*norm(dk)^2
    % model=2 means -F(zk)'*dk ≥ sigma*tk*norm(F(zk))*norm(dk)^2
    % model=3 means -F(zk)'*dk ≥ sigma*tk*norm(F(zk))/(1+norm(F(zk)))*norm(dk)^2
    % model=4 means -F(zk)'*dk ≥ sigma*tk*max(lambda,min(nu,norm(Fz_new,2)))*norm(dk)^2
    if model==1
        t = gamma;
        z0_new = v0+t*dk;
        Fz0_new = feval(nprob,z0_new,1);
        NF = NF+1;
        Fz0_newtdk = -Fz0_new'*dk;
        % check the Armijo-type line search condition
        while Fz0_newtdk < sigma*t*Normdk2 %&& t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            z0_new = v0+t*dk;
            Fz0_new = feval(nprob,z0_new,1);
            NF = NF+1;
            Fz0_newtdk = -Fz0_new'*dk;
        end %%% End Armijo-type line search %%%
        NormFz0_new = norm(Fz0_new);
    elseif model==2
        t = gamma;
        z0_new = v0+t*dk;
        Fz0_new = feval(nprob,z0_new,1);
        NF = NF+1;
        NormFz0_new = norm(Fz0_new);
        Fz0_newtdk = -Fz0_new'*dk;
        % check the Armijo-type line search condition
        while Fz0_newtdk < sigma*t*NormFz0_new*Normdk2 %&& t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            z0_new = v0+t*dk;
            Fz0_new = feval(nprob,z0_new,1);
            NF = NF+1;
            NormFz0_new = norm(Fz0_new);
            Fz0_newtdk = -Fz0_new'*dk;
        end %%% End Armijo-type line search %%%
    elseif model==3
        t = gamma;
        z0_new = v0+t*dk;
        Fz0_new = feval(nprob,z0_new,1);
        NF = NF+1;
        NormFz0_new = norm(Fz0_new);
        Fz0_newtdk = -Fz0_new'*dk;
        % check the Armijo-type line search condition
        while Fz0_newtdk < sigma*t*NormFz0_new/(1+NormFz0_new)*Normdk2 && t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            z0_new = v0+t*dk;
            Fz0_new = feval(nprob,z0_new,1);
            NF = NF+1;
            NormFz0_new = norm(Fz0_new);
            Fz0_newtdk = -Fz0_new'*dk;
        end %%% End Armijo-type line search %%%
    else
        t = gamma;
        z0_new = v0+t*dk;
        Fz0_new = feval(nprob,z0_new,1);
        NF = NF+1;
        NormFz0_new = norm(Fz0_new);
        Fz0_newtdk = -Fz0_new'*dk;
        % check the Armijo-type line search condition
        while Fz0_newtdk < sigma*t*min(0.5,0.2*NormFz0_new)*Normdk2 %&& t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            z0_new = v0+t*dk;
            Fz0_new = feval(nprob,z0_new,1);
            NF = NF+1;
            NormFz0_new = norm(Fz0_new);
            Fz0_newtdk = -Fz0_new'*dk;
        end %%% End Armijo-type line search %%%
    end 
    Fz0 = Fz0_new;
    NormFz0 = NormFz0_new;
%     if NormFz0<=epsilon
%         L1 = 1;
%         NormF = NormFz0; % the final norm of equations
%         break;
%     end
    xik = t*Fz0_newtdk/NormFz0^2;
    z1 = v0-rho*xik*Fz0;
    % compute the next iteration 
    x1 =z1;% feval(nprob,z1,2);
    Fx1 = feval(nprob,x1,1);
    NF = NF+1;
    NormFx1 = norm(Fx1);
    if NormFx1<=epsilon
        L1 = 1;
        NormF = NormFx1;
        break;
    end
    
    % update the iteration
    d_k_prev =dk;
    x0_old = x0;
    x0 = x1;
    y0_old=v0;
    NormFx0 = NormFx1;
end
if L1==1
    Itr = k;
    Tcpu = toc;
else
    NF = NaN;
    Itr = NaN;
    Tcpu = NaN;
    NormF = NaN;
end
