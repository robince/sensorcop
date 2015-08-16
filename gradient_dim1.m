function d = gradient_dim1(f)
% take central difference gradient across the first dimension
% no spacing

d = zeros(size(f),class(f));
n = size(f,1);

% Take forward differences on left and right edges
d(1,:) = f(2,:) - f(1,:);
d(n,:) = f(n,:) - f(n-1,:);

% Take centered differences on interior points
d(2:n-1,:) = (f(3:n,:)-f(1:n-2,:));