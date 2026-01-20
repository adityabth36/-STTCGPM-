
function out= problem3(x,mode)
%Problem 2 in ¡¶Spectral Gradient Projection Method for solving Nonlinear
% Equations with Convex Constraints¡·, Z.S. Yu, J. Lin, J. Sun, Y.H. Xiao, L.Y. Liu, Z.H. Li 2009

n = length(x);
if mode==1               % compute F(x)
  Fx=ones(n,1);
  for i=1:n
     Fx(i)=min(min(abs(x(i)+1), x(i)), max(abs(x(i)-1), x(i)^3));
  end
  out=Fx; 
elseif  mode==2          % compute the projection
     out=max(x,0);
end
end
    
        
    
    


