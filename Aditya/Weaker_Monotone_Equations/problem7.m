
function out= problem7(x)

% n=length(x);
% if mode==1
%     a=5;
%   Fx=ones(n,1);
%   for i=1:n-1
%     Fx(i)=(x(i)-1)*(1 + a*(x(i+1)-1)^2) ;
%   end
%   Fx(n)=(x(n)-1)*(1+a*(x(1)-1)^2); 
%   out=Fx;
% elseif  mode==2
%     out=max(x,0);
% end
% end
    
    a=5;
  Fx=ones(2,1);
    Fx(1)=(x(1)-1)*(1 + a*(x(2)-1)^2) ;
    Fx(2)=(x(2)-1)*(1+a*(x(1)-1)^2); 
  out=Fx;
end        
    
    


