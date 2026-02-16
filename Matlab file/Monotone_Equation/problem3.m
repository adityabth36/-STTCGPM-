
function out= problem3(x,mode)

n = length(x);
if mode==1
  Fx=ones(n,1);
    % Fx(1)=exp(x(1))-1;
  for i=1:n-1
    Fx(i)=exp(x(i)) + x(i)-1 ;
  end
      Fx(n)=2*x(n)-x(n-1)+exp(x(n))-1;
  out=Fx;
elseif  mode==2
    out=max(x,0);
end

end
    
        
    
    


