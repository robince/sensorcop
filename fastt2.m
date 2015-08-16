function t = fastt2(X, stim)

X1 = X(stim==0,:);
X2 = X(stim==1,:);

n1 = size(X1,1);
n2 = size(X2,1);

mX1 = mean(X1,1);
mX2 = mean(X2,1);

vX1 = var(X1,0,1);
vX2 = var(X2,0,1);

% tnum = mean(X1,1) - mean(X2,1);
% % equal variances
% S12 = sqrt(((n1-1)*vX1 + (n2-1)*vX2) ./ (n1+n2-2));
% t = (mX1 - mX2) ./ (S12*sqrt( (1/n1) + (1/n2) ));
% unequal variances
S12 = sqrt(vX1/n1 + vX2/n2);
t = (mX1 - mX2) ./ S12;