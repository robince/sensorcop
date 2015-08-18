function t2 = fastht22(X, stim)
% two sample Hotelling T2 test with unequal covariance

X1 = X(stim==0,:,:);
X2 = X(stim==1,:,:);

n1 = size(X1,1);
n2 = size(X2,1);

mX1 = mean(X1,1);
mX2 = mean(X2,1);

% demean
cX1 = bsxfun(@minus, X1, mX1);
cX2 = bsxfun(@minus, X2, mX2);

% cov with mtimesx

% covX1 = mtimesx(cX1,'T',cX1,'N');
% covX2 = mtimesx(cX2,'T',cX2,'N');

covX1 = mmx('square',cX1,[],'t') ./ (n1-1);
covX2 = mmx('square',cX2,[],'t') ./ (n2-1);

% pool covariance estimates (homogenous)
% Cp = (covX1 + covX2) / (n1 + n2 - 2);
% separate cov (heterogenous)
Cp = (covX1./n1) + (covX2./n2);
dm = mX1-mX2;
dm = reshape(dm, [2 1 size(dm,3)]);

invCpdm = mmx('backslash', Cp, dm);
t2 = squeeze(mmx('mult',dm,invCpdm,'tn'));
