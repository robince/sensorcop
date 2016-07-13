function Zb = nonoverlappingBB(Z,nboot,sim,L)
% PURPOSE: Stationary Block Bootstrap for a vector time series
% ------------------------------------------------------------
% SYNTAX: Zb = stationaryBB(Z,sim,L);
% ------------------------------------------------------------
% OUTPUT: Zb : nxkz resampled time series
% ------------------------------------------------------------
% INPUT:  Z   : nxkz --> vector time series to be resampled
%         sim : 1x1  --> type of bootstrap: 
%                       1 => stationary geometric pdf
%                       2 => stationary uniform pdf
%                       3 => circular (non-random b)
%         L --> block size, depends on sim
%           sim = 1 --> L:1x1 expected block size
%           sim = 2 --> L:2x1 lower and upper limits for uniform pdf
%           sim = 3 --> L:1x1 fixed block size (non-random b)
%        If L=1 and sim=3, the standard iid bootstrap is applied
% ------------------------------------------------------------
% LIBRARY: loopBB [internal]
% ------------------------------------------------------------
% SEE ALSO: overlappingBB, seasBB
% ------------------------------------------------------------
% REFERENCES: Politis, D. and Romano, J. (1994) "The starionary 
% bootstrap", Journal of the American Statistical Association, vol. 89,
% n. 428, p. 1303-1313.
% Politis, D. and White, H (2003) "Automatic block-length
% selection for the dependent bootstrap", Dept. of Mathematics, Univ.
% of California, San Diego, Working Paper.
% ------------------------------------------------------------

% written by:
%  Enrique M. Quilis
%  Macroeconomic Research Department
%  Fiscal Authority for Fiscal Responsibility (AIReF)
%  <enrique.quilis@airef.es>

% Version 2.1 [October 2015]

% Dimension of time series to be bootstrapped
% n specified by user
[n,kz] = size(Z);

% ------------------------------------------------------------
%  ALLOCATION
% ------------------------------------------------------------
I = zeros(1,nboot);
b = zeros(1,nboot);
xb = -999.99*ones(nboot,1);

% ------------------------------------------------------------
% WRAPPING THE TIME SERIES AROUND THE CIRCLE
% ------------------------------------------------------------
Z = [Z; Z(1:n-1,:)];

% ------------------------------------------------------------
% INDEX SELECTION
% ------------------------------------------------------------
% I = round(1+(n-1)*rand(1,nboot));

% ------------------------------------------------------------
% BLOCK SELECTION
% ------------------------------------------------------------
% switch sim
% case 1 % Stationary BB, geometric pdf
%    b = geornd(1/L(1),1,nboot);
% case 2 % Stationary BB, uniform pdf   
%    b = round(L(1)+(L(2)-1)*rand(1,nboot));
% case 3 % Circular bootstrap (fixed block size)
%    b = L(1) * ones(1,nboot);
% end

nblock = ceil(n ./ L);
block_start = ((0:(nblock-1))*L)+1;
block_order = randperm(nblock);
I = block_start(block_order);
b = L(1) * ones(1,nboot);

% ------------------------------------------------------------
% BOOTSTRAP REPLICATION
% ------------------------------------------------------------
Zb = [];
for j=1:kz
   Zb = [Zb loopBB(Z(:,j),nboot,b,I)];
end
if size(Zb,1) ~= nboot
    error('not enough samples, try smaller block size')
end

% ============================================================
% loopBB ==> UNIVARIATE BOOTSTRAP LOOP
% ============================================================
function xb = loopBB(x,n,b,I);

h=1;
% xb = zeros(n,1);
for m=1:length(I)
   for j=1:b(m)
      xb(h) = x(I(m)+j-1);
      h = h + 1; 
      if (h == n+1); break; end;
   end
   if (h == n+1); break; end;
end

xb=xb';