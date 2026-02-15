
function out= problem7(x,mode)

n=length(x);
if mode==1
    a=5;
  Fx=ones(n,1);
  for i=1:n-1
    Fx(i)=(x(i)-1)*(1 + a*(x(i+1)-1)^2) ;
  end
  Fx(n)=(x(n)-1)*(1+a*(x(1)-1)^2); 
  out=Fx;
elseif  mode==2
    r=2*sqrt(n);
    nrm = vecnorm(x);       % Euclidean norm of each column
    scale = min(1, r ./ nrm);
    out = x .* scale;         % columnwise scaling
    % out=max(x,0);
end
% end
% elseif mode ==2   
% out = min(max(x, -2), 2);
% end
end        

    
