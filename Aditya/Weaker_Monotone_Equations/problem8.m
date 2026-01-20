function out= problem8(x,mode)
%Example 1 in ¡¶Solving pseudomonotone variational inequalities and
% pseudoconvex optimization problems using the projection neural network¡·,Hu, X., Wang, J.,
% IEEE Trans. Neural Networks 2006, 17(6): 1487-1499

n=length(x);
% if mode==1
%     v = x - ones(n,1);  % shift by (1,...,1)
%     % --- Construct D(x) diagonal matrix ---
%     d = zeros(n,1);
%     for i = 1:n
%         ip1 = mod(i,n) + 1; % cyclic index i+1
%         d(i) = 1 + exp(-(x(i)-1)^2) ;%+ log(1 + (x(ip1)-1)^2);
%     end
%     D = diag(d);
%     % --- Construct B(x) skew-symmetric matrix ---
%     B = zeros(n);
%     for i = 1:n
%         ip1 = mod(i,n)+1 ; % cyclic index i+1
%         B(i,ip1) = sin(x(i) - x(ip1));
%         B(ip1,i) = -B(i,ip1);
%     end
%     Fx = (D + B) * v;
%     out=Fx;
% elseif  mode==2
%     out=max(x,0);
% end
% end

%%  Optimize way to construct diagonal and skew symmetric matrix
if mode==1
     v = x - 1;  % shift by (1,...,1)
    % --- Diagonal part D*v ---
    ip1 = [2:n 1];  
    im1 = [n 1:n-1];
    d = 1 + exp(-(x-1).^2)+ log(1 + (x(ip1)-1).^2);
    Dv = d .* v;
    % --- Skew-symmetric part B*v ---
    % ip1 = [2:n 1];                        % cyclic index shift
    s = sin(x - x(ip1));                  % B(i,ip1)
    Bv = s .* v(ip1) - s(im1) .* v(im1); 
    % --- Final result ---
    out = Dv + Bv;
elseif  mode==2
    out=max(x,0);
end
end



