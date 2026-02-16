% init 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [nprob,n,m,x0]=initf(nprob)
% This function sets n,m, and the standard starting    
% point based on the nprob and returns it to initpt     
% function.                                                                                                    
% Created on 10/30/94 by Madhu Lamba                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nprob,n] = init(NO)

%global FIRSTIME

switch NO
   %% nprob='problem2';
      case 1
        nprob='problem9';
        n=1000;
    
    %% nprob='problem2';
      case 2
        nprob='problem9';
        n=3000;
        
    %% nprob='problem2';
      case 3
        nprob='problem9';
        n=5000;
        
    %% nprob='problem2';
      case 4
        nprob='problem9';
        n=8000;
    
    %% nprob='problem2';
      case 5
        nprob='problem9';
        n=10000;
        
    %% nprob='problem2';
      case 6
        nprob='problem9';
        n=30000;
    
    %% nprob='problem2';
      case 7
        nprob='problem9';
        n=50000;
        
    %% nprob='problem2';
      case 8
        nprob='problem9';
        n=80000;
        
    %% nprob='problem2';
      case 9
        nprob='problem9';
        n=100000;     
end