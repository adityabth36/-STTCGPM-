function out= problem9(x,mode)
%Example 4.1 in ¡¶A Projection Method for a system of Nonlinear Monotone
%Equations with convex constraints¡·, C. Wang, Y. Wang, C. Xu, 2007

n = length(x);
if mode==1               % compute F(x)
  Fx=ones(n,1);
  for i=1:n
   Fx(i)=log(abs(x(i)+1))-x(i)/n ;
    
  end
  out=Fx;
elseif  mode==2          % compute the projection
    out=max(x,0);
end