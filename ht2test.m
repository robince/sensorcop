function t2 = ht2test(x,y)
% two sample Hotelling T2 test with equal covariance
mx = mean(x,1);
my = mean(y,1);

Nx = size(x,1);
Ny = size(y,1);

xn = bsxfun(@minus, x, mx);
yn = bsxfun(@minus, y, my);
Cx = xn'*xn / (Nx-1);
Cy = yn'*yn / (Ny-1);


% pool covariance estimates (homogenous)
Cp = ((Nx-1)*Cx + (Ny-1)*Cy) / (Nx + Ny - 2);
% separate cov (heterogenous)
% Cp = (Cx./Nx) + (Cy./Ny);
dm = mx-my;
t2 = dm * inv(Cp) * dm';
