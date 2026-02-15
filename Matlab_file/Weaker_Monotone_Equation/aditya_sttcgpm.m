
% This is a demo of DFPM for solving constrained nonlinear weaker-monotone equations
% equations of the form
%   F(x)=0, x\in C, 
% where C is a nonempty closed convex set.

clc;
clear all
close all

% set random number seed
rng(2025)

ITR_max = 2000;
% setup TXT document
fid_tex=fopen('mytext.txt','w'); 
% problem_set = 1:72;
 problem_set =1:27;
% set parameters
np = length(problem_set); % from problem 1 to problem 8 % the number of the test problems
ns = 1;   % the number of the test algorithms
T = zeros(np,ns);
F = zeros(np,ns);
N = zeros(np,ns);
G = zeros(np,ns);

% set parameters for inertial derivative-free projection methods (IDFPMs)
% set parameters for STTCGPM
para1.Itr_max = ITR_max;
para1.gamma = 0.8;         % the initial guess
para1.sigma =  0.01;%10^-4;         % the coefficient of line search 
para1.tau = 0.6;         % the compression ratio
para1.alpha = 0.2;       % the coefficient of inertial step
para1.rho = 1.9;         % the relaxation factor 


% set parameters for inertial projection methods (FITTCGPM-PRP/FITTCGPM-DY/FITTCGPM-our)
para2.Itr_max = ITR_max;
para2.gamma = 0.8;         % the initial guess
para2.sigma = 0.01;         % the coefficient of line search 
para2.tau = 0.6;         % the compression ratio
para2.alpha = 0.2;       % the coefficient of inertial step
para2.rho = 1.8;         % the relaxation factor 

% set parameters for GITDFPA
para3.Itr_max = ITR_max;
para3.gamma = 0.8;         % the initial guess
para3.sigma = 0.01;         % the coefficient of line search 
para3.tau = 0.6;         % the compression ratio
para3.alpha = 0.2;       % the coefficient of inertial step
para3.rho = 1;         % the relaxation factor 

% run
for index=1:9  %np
    Num = problem_set(index);
    [name,n] = init(Num);
    progress_r = [];
    for repeats = 1:5
             x0 = 5*rand(n,1);    
     
         [T1,NFF1,NI1,G1] = DFPM(Num,'STTCGPM',1,2,x0,para1); % acceleration
        progress_r = [progress_r;NI1,NFF1,T1,G1];
    end
    TM = mean(progress_r);
    fprintf(fid_tex,'%s %d & %.1f/%.1f/%.3f/%.2e & %.1f/%.1f/%.3f/%.2e\n& %.1f/%.1f/%.3f/%.2e & %.1f/%.1f/%.3f/%.2e\\\\ \r\n', ... 
                name,n,TM);
    T(index,:) = [TM(3)];
    F(index,:) = [TM(2)];
    N(index,:) = [TM(1)];
    G(index,:) = [TM(4)];
end
