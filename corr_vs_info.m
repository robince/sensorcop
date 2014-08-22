% corrleation vs information for two gaussian variables
C = linspace(-1,1,100);
NC = length(C);

I = zeros(size(C));
for i=1:NC
    
    thsI = log(2*pi*exp(1)) - 1 - log(2*pi) - 0.5*log(det([ 1 C(i); C(i) 1]));
    I(i) = thsI;
end

plot(C,I)
xlabel('Correlation')
ylabel('Information (bits)')
title('Two Gaussian Variables')