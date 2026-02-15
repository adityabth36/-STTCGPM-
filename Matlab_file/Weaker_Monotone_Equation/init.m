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
    %% nprob='problem1';
      case 1
        nprob='problem7';
        n=1000;
    
    %% nprob='problem1';
      case 2
        nprob='problem7';
        n=3000;
        
    %% nprob='problem1';
      case 3
        nprob='problem7';
        n=5000;
        
    %% nprob='problem1';
      case 4
        nprob='problem7';
        n=8000;
    
    %% nprob='problem1';
      case 5
        nprob='problem7';
        n=10000;

    %% nprob='problem1';
      case 6
        nprob='problem7';
        n=30000;
    
    %% nprob='problem1';
      case 7
        nprob='problem7';
        n=50000;
        
    %% nprob='problem1';
      case 8
        nprob='problem7';
        n=80000;
        
    %% nprob='problem1';
      case 9
        nprob='problem7';
        n=100000;
        
    %% nprob='problem2';
      case 10
        nprob='problem8';
        n=1000;
    
    %% nprob='problem2';
      case 11
        nprob='problem8';
        n=3000;
        
    %% nprob='problem2';
      case 12
        nprob='problem8';
        n=5000;
        
    %% nprob='problem2';
      case 13
        nprob='problem8';
        n=8000;
    
    %% nprob='problem2';
      case 14
        nprob='problem8';
        n=10000;
        
    %% nprob='problem2';
      case 15
        nprob='problem8';
        n=30000;
    
    %% nprob='problem2';
      case 16
        nprob='problem8';
        n=50000;
        
    %% nprob='problem2';
      case 17
        nprob='problem8';
        n=80000;
        
    %% nprob='problem2';
      case 18
        nprob='problem8';
        n=100000;
   %% nprob='problem2';
      case 19
        nprob='problem9';
        n=1000;
    
    %% nprob='problem2';
      case 20
        nprob='problem9';
        n=3000;
        
    %% nprob='problem2';
      case 21
        nprob='problem9';
        n=5000;
        
    %% nprob='problem2';
      case 22
        nprob='problem9';
        n=8000;
    
    %% nprob='problem2';
      case 23
        nprob='problem9';
        n=10000;
        
    %% nprob='problem2';
      case 24
        nprob='problem9';
        n=30000;
    
    %% nprob='problem2';
      case 25
        nprob='problem9';
        n=50000;
        
    %% nprob='problem2';
      case 26
        nprob='problem9';
        n=80000;
        
    %% nprob='problem2';
      case 27
        nprob='problem9';
        n=100000;     
end