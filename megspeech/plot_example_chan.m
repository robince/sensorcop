load(fullfile(data_dir,'sub2_dat.mat'));
Fs = 50;

[Npln Nt] = size(plndat);

% pln sum
Nch = Npln / 2;
if round(Nch) ~= Nch
    error('Odd number of planar channels')
end
% planar sum square
meg = sqrt(plndat(1:Nch,:).^2 + plndat((Nch+1):end,:).^2);
meg = meg';

%%
megH = plndat(1:Nch,:)';
megV = plndat((Nch+1):end,:)';
megdH = gradient_dim1(megH);
megdV = gradient_dim1(megV);
meg2d = permute(cat(3,megH,megV),[1 3 2]);
meg4d = permute(cat(3,megH,megV,megdH,megdV),[1 3 2]);

%%

chi = 231;
thsNblock = 45;

delays = 0:40;
Ndel = length(delays);

Nthread = 4;
blocklen = 10 * Fs; % 10 s blocks
Nblock = ceil(Nt / blocklen);


Icop = zeros(1,Ndel);
Ib2 = zeros(1,Ndel);
Ib4 = zeros(1,Ndel);
Ib8 = zeros(1,Ndel);
sp = zeros(1,Ndel);
pe = zeros(1,Ndel);


blkidx = randperm(Nblock, thsNblock);

% quantise blocks used in this rep
bmeg = cell(1,thsNblock);
bspc = cell(1,thsNblock);
blen = zeros(thsNblock,1);
bmeg2 = cell(1,thsNblock);
bmeg4 = cell(1,thsNblock);
for bi=1:thsNblock
    idx = block_index(blkidx(bi),blocklen,Nt);
    bmeg{bi} = meg(idx,chi);
    bspc{bi} = fltspc(idx);
    blen(bi) = length(idx);
    bmeg2{bi} = meg2d(idx,:,chi);
    bmeg4{bi} = meg4d(idx,:,chi);
end

thsmeg = cell2mat(bmeg');
thsmeg2 = cell2mat(bmeg2');
[c s] = pca(thsmeg2);
thspca = s(:,1);
thsmeg4 = cell2mat(bmeg4');
thsspc = cell2mat(bspc');

cmeg = copnorm_slice_omp_c_double(thsmeg, Nthread);
cpca = copnorm_slice_omp_c_double(thspca, Nthread);
cmeg2 = copnorm_slice_omp_c_double(thsmeg2, Nthread);
cmeg4 = copnorm_slice_omp_c_double(thsmeg4, Nthread);

qmeg2 = int16(bin.eqpop(thsmeg,2));
qmeg4 = int16(bin.eqpop(thsmeg,4));
qmeg8 = int16(bin.eqpop(thsmeg,8));

cspc = copnorm(thsspc);
qspc2 = int16(bin.eqpop(thsspc,2));
qspc4 = int16(bin.eqpop(thsspc,4));
qspc8 = int16(bin.eqpop(thsspc,8));

bcmeg = cell(1,thsNblock);
bcmeg2 = cell(1,thsNblock);
bcmeg4 = cell(1,thsNblock);
bcpca = cell(1,thsNblock);
bqmeg2 = cell(1,thsNblock);
bqmeg4 = cell(1,thsNblock);
bqmeg8 = cell(1,thsNblock);
bcspc = cell(1,thsNblock);
bqspc2 = cell(1,thsNblock);
bqspc4 = cell(1,thsNblock);
bqspc8 = cell(1,thsNblock);
ii = 1;
for bi=1:thsNblock
    idx = ii:(ii+blen(bi)-1);
    ii = ii+blen(bi);
    bcmeg{bi} = cmeg(idx,:);
    bcmeg2{bi} = cmeg2(idx,:);
    bcmeg4{bi} = cmeg4(idx,:);
    bcpca{bi} = cpca(idx,:);
    
    bqmeg2{bi} = qmeg2(idx,:);
    bqmeg4{bi} = qmeg4(idx,:);
    bqmeg8{bi} = qmeg8(idx,:);
    
    bcspc{bi} = cspc(idx);
    bqspc2{bi} = qspc2(idx);
    bqspc4{bi} = qspc4(idx);
    bqspc8{bi} = qspc8(idx);
end

% TRUE VALUE
for di=1:Ndel
    % align
    d = delays(di);
    
    [dmeg, dspc] = block_delay(bcmeg, bcspc, d);
    Icop(di) = info_cc_slice_nobc_omp(reshape(dmeg,[size(dmeg,1) 1 size(dmeg,2)]),1, dspc, size(dmeg,1),Nthread);
    
    [dqmeg, dqspc] = block_delay(bqmeg2, bqspc2, d);
    Ib2(di) = info.calc_info_slice_omp_integer_c_int16_t(dqmeg,2,dqspc,2,size(dqmeg,1),Nthread);
    [dqmeg, dqspc] = block_delay(bqmeg4, bqspc4, d);
    Ib4(di) = info.calc_info_slice_omp_integer_c_int16_t(dqmeg,4,dqspc,4,size(dqmeg,1),Nthread);
    [dqmeg, dqspc] = block_delay(bqmeg8, bqspc8, d);
    Ib8(di) = info.calc_info_slice_omp_integer_c_int16_t(dqmeg,8,dqspc,8,size(dqmeg,1),Nthread);
    
    Ib2(di) = Ib2(di) - 1./(2*size(dqmeg,1)*log(2));
    Ib4(di) = Ib4(di) - 9./(2*size(dqmeg,1)*log(2));
    Ib8(di) = Ib8(di) - 49./(2*size(dqmeg,1)*log(2));
    
    
    [dmeg, dspc] = block_delay(bmeg, bspc, d);
    sp(di) = (corr(dspc, dmeg, 'type', 'spearman'));
    pe(di) = (corr(dspc, dmeg, 'type', 'pearson'));
    %             fres.ke(:,di) = corr(dspc, dmeg, 'type', 'kendall');
end

%%
figure
n = 3;
dx = delays*20;

subplot(n,1,1)
plot(dx,Icop)
axis tight

subplot(n,1,2)
plot(dx,Ib2)
hold all
plot(dx,Ib4)
plot(dx,Ib8)
axis tight

subplot(n,1,3)
plot(dx,sp)
axis tight
xlabel('Lag (ms)')

set(gcf,'Pos',[680   859   286   239])

% subplot(414)
% plot(pe)

% figure
% plot(pe)

%%
% TRUE VALUE
for di=1:Ndel
    % align
    d = delays(di);
    
    [dmeg, dspc] = block_delay(bcmeg2, bcspc, d);
    Icop2d(di) = info_gg(dmeg,dspc,false,false,false);
    IcopH(di) = info_gg(dmeg(:,1),dspc,false,false,false);
    IcopV(di) = info_gg(dmeg(:,2),dspc,false,false,false);
    spH(di) = (corr(dspc, dmeg(:,1), 'type', 'spearman'));
    spV(di) = (corr(dspc, dmeg(:,2), 'type', 'spearman'));
    [dmeg, dspc] = block_delay(bcmeg4, bcspc, d);
    Icop4d(di) = info_gg(dmeg,dspc,false,false,false);
    [dmeg, dspc] = block_delay(bcpca, bcspc, d);
    Icoppc1(di) = info_gg(dmeg,dspc,false,false,false);
end

%%
figure
plot(dx,Icop)
hold all
% figure
plot(dx,Icop2d)
plot(dx,Icop4d)
%%
figure
subplot(2,1,1)
plot(dx,Icop)
hold all
plot(dx,Icop2d,'g')


subplot(2,1,2)
plot(dx,spH)
hold all
plot(dx,spV)

% plot(Icoppc1,'k')


%%
figure
subplot(211)
plot(IcopH)
hold all
plot(IcopV)
subplot(212)
plot(spH)
hold all
plot(spV)


