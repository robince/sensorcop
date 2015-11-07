% corrleation vs information for two gaussian variables
C = linspace(-1,1,100);
NC = length(C);

I = zeros(size(C));
for i=1:NC
    
    thsI = log(2*pi*exp(1)) - 1 - log(2*pi) - 0.5*log(det([ 1 C(i); C(i) 1]));
    I(i) = thsI;
end


%%
figure
plot(C,I,'k','LineWidth',2)
xlabel('Correlation','FontSize',16)
ylabel('Information (bits)','FontSize',16)
title('Two Gaussian Variables','FontSize',16)
xlim([-1.1 1.1])

vline(-1,'k:')
vline(1,'k:')
axis square
set(gca,'FontSize',16)
set(gca,'TickLength',[0.02 0.025])
set(gca,'TickDir','out')
box off
