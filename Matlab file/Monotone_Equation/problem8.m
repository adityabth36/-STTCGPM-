function out= problem8(x,mode)
n = length(x);
if mode==1               % compute F(x)
  Fx=ones(n,1);
  for i=1:n
   Fx(i)=2*x(i)-sin(abs(x(i)));
    
  end
  out=Fx;
elseif  mode==2          % compute the projection
    out=max(x,0);
end
