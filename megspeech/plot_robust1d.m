% plot results from simulation

% samp_size = [25 50 100 200 500];
% stats = {'Icop' 'Ib2' 'Ib4' 'Ib8' 't' 'ks'};

Ncorr = length(corrupt_prct);
% stats = {'Icop' 'Ib2' 'Ib4' 'Ib8' 't' 'ks'};
Nstats = length(stats);
samp_size = corrupt_prct;
%%
Nrep = size(res,1);

cm = [];
for si=1:Nstats
    stat = stats{si};
    t = cellfun(@(x) x.(stat).cm, res, 'Uniform',false);
    cm.(stat) = cell2mat(reshape(t,[1 1 Nrep Ncorr]));
end


%% SENSITIVITY
figure
% hold all

for si=1:Nstats
    st = cm.(stats{si});
    % sensivity
    tp = squeeze(st(1,1,:,:));
    tn = squeeze(st(2,2,:,:));
    % TODO: check
    fp = squeeze(st(1,2,:,:));
    fn = squeeze(st(2,1,:,:));
    
    sens = tp ./ (tp+fn);
    
    plot(10*samp_size, mean(sens,1));
%     semilogx(samp_size, mean(sens,1));
    hold all
end
legend(stats)
title('Sensitivity')
ylim([0 0.5])
xlabel('Recording (s)')

%% PRECISION
figure
hold all

for si=1:Nstats  
    st = cm.(stats{si});
    % sensivity
    tp = squeeze(st(1,1,:,:));
    tn = squeeze(st(2,2,:,:));
    % TODO: check
    fp = squeeze(st(1,2,:,:));
    fn = squeeze(st(2,1,:,:));
    
%     sens = tp ./ (tp+fn);
    precision = tp ./ (tp + fp);
    
    plot(samp_size, nanmean(precision,1));
end
legend(stats)
title('Precision')
% ylim([0.9 1])

%% SPECIFICITY
figure
hold all

for si=1:Nstats  
    st = cm.(stats{si});
    % sensivity
    tp = squeeze(st(1,1,:,:));
    tn = squeeze(st(2,2,:,:));
    % TODO: check
    fp = squeeze(st(1,2,:,:));
    fn = squeeze(st(2,1,:,:));
    
%     sens = tp ./ (tp+fn);
%     precision = tp ./ (tp + fp);
    specificity = tn ./ (tn + fp);
    
    plot(samp_size, mean(specificity,1));
end
legend(stats)
title('Specificity')
ylim([0.95 1.05])
box on


%% MATTHEWS
figure
% hold all

for si=1:Nstats
    st = cm.(stats{si});
    % sensivity
    tp = squeeze(st(1,1,:,:));
    tn = squeeze(st(2,2,:,:));
    % TODO: check
    fp = squeeze(st(1,2,:,:));
    fn = squeeze(st(2,1,:,:));
    
    denom = sqrt((tp+fp).*(tp+fn).*(tn+fp).*(tn+fn));
    mcc = (tp.*tn - fp.*fn) ./ denom;
    
    plot(samp_size, nanmean(mcc,1));
%     semilogx(samp_size, nanmean(mcc,1));
    hold all
end
legend(stats,'Location','SouthEast')
title('MCC')
% ylim([0 1])

%% INFORMATION
figure


ent = @(p) sum(-p(p>0).*log2(p(p>0)));

for si=1:Nstats
    st = cm.(stats{si});
    % sensivity
    tp = squeeze(st(1,1,:,:));
    tn = squeeze(st(2,2,:,:));
    % TODO: check
    fp = squeeze(st(1,2,:,:));
    fn = squeeze(st(2,1,:,:));
    
    I = zeros(size(tp));
    for ci=1:size(tp(:))
        Pj = zeros(4,1);
        N = tp(ci) + fp(ci) + fn(ci) + tn(ci);
        Pj(1) = tp(ci);
        Pj(2) = fp(ci);
        Pj(3) = tn(ci);
        Pj(4) = fn(ci);
        Pj = Pj ./ N;
        
        Hj = ent(Pj);
        
        Ptrue = zeros(2,1);
        Ptrue(1) = tp(ci) + fn(ci);
        Ptrue(2) = fp(ci) + tn(ci);
        Ptrue = Ptrue ./ N;
        Htrue = ent(Ptrue);
        
        Ptest = zeros(2,1);
        Ptest(1) = tp(ci) + fp(ci);
        Ptest(2) = fn(ci) + tn(ci);
        Ptest = Ptest ./ N;
        Htest = ent(Ptest);
        I(ci) = (Htrue + Htest - Hj);
        [i j] = ind2sub(size(tp),ci);
        I(ci) = (Htrue + Htest - Hj) - 1./(2*N*log(2));
        I(ci) = I(ci) ./ max(Htrue,Htest);
    end
    
%     mean(I,1)
    plot(samp_size, mean(I,1));
%     semilogx(samp_size, mean(I,1));
    hold all
end
legend(stats,'Location','SouthEast')
title('MI')