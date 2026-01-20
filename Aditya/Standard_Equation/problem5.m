function out= problem5(x,mode)

 n = length(x);
if mode==1               % compute F(x)
  Fx=ones(n,1);
    Fx(1) = 2*x(1) ;
  for i=2:n
      Fx(i) = 0.5*cos(x(i-1))+2*x(i)-2 ;
  end
  out=Fx;
elseif  mode==2          % compute the projection
    out=max(x,0);
end  
end
    
        
    
    


