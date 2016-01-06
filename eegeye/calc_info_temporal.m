% calculate temporal info for RE/LE right/left eye

% flti = 1; % regular non-causal
flti = 2; % causal
Nele = 2;
Nstm = 2;

%%
% subi = 1;
parfor subi=1:Ns
    subid = subjects{subi};
    fname = sprintf('data_%s_%s.mat',subid,fqnames{flti});
    dat = matfile(fullfile(data_dir,fname));
    eyestm = dat.eyebubs;
    otl = dat.ferp(:,:,dat.LEmaxMI);
    otr = dat.ferp(:,:,dat.REmaxMI);
    fdat = cat(3, otl, otr);
    
    ddat = zeros(size(fdat));
    d2dat = zeros(size(fdat));
    for eli=1:2
        for trli=1:size(fdat,1)
            ddat(trli,:,eli) = gradient(squeeze(fdat(trli,:,eli)),2);
            d2dat(trli,:,eli) = gradient(ddat(trli,:,eli),2);
        end
    end
    
    cstm = copnorm(eyestm);
    cdat = copnorm(fdat);
    cddat = copnorm(ddat);
    cd2dat = copnorm(d2dat);
    
    [Ntrl, Nt, ~] = size(fdat);
    res = [];
    
    %
    % TIMECOURSES
    %
    
    % conditional info timecourses
    res.p.Itc = zeros(Nt,2,2);
    res.pv.Itc = zeros(Nt,2,2);
    res.pva.Itc = zeros(Nt,2,2);
    % direct info timecourses
    res.p.It = zeros(Nt,2,2);
    res.pv.It = zeros(Nt,2,2);
    res.pva.It = zeros(Nt,2,2);
    for ti=1:Nt
        for stmi=1:2
            for eleci=1:2
                otherstm = 1 + mod(stmi,2);
                
                res.pva.Itc(ti,stmi,eleci) = cmi_ggg(cstm(:,stmi), [cdat(:,ti,eleci) cddat(:,ti,eleci) cd2dat(:,ti,eleci) ], cstm(:,otherstm), true, false, false);
                res.pva.It(ti,stmi,eleci) = info_gg(cstm(:,stmi), [cdat(:,ti,eleci) cddat(:,ti,eleci) cd2dat(:,ti,eleci) ], true, false, false);
                
                res.pv.Itc(ti,stmi,eleci) = cmi_ggg(cstm(:,stmi), [cdat(:,ti,eleci) cddat(:,ti,eleci) ], cstm(:,otherstm), true, false, false);
                res.pv.It(ti,stmi,eleci) = info_gg(cstm(:,stmi), [cdat(:,ti,eleci) cddat(:,ti,eleci) ], true, false, false);
                
                res.p.Itc(ti,stmi,eleci) = cmi_ggg(cstm(:,stmi), [cdat(:,ti,eleci) ], cstm(:,otherstm), true, false, false);
                res.p.It(ti,stmi,eleci) = info_gg(cstm(:,stmi), [cdat(:,ti,eleci) ], true, false, false);
                
            end
        end
    end
    
    %
    % SELF TEMPORAL INTERACTION
    %
    
    xl = [0 400];
    idx = (time>xl(1)) & (time<xl(2));
    intNt = sum(idx);
    inttime = time(idx);
    % conditional
    res.p.intItc = zeros(intNt,intNt,2,2);
    res.pv.intItc = zeros(intNt,intNt,2,2);
    res.pva.intItc = zeros(intNt,intNt,2,2);
    % direct
    res.p.intIt = zeros(intNt,intNt,2,2);
    res.pv.intIt = zeros(intNt,intNt,2,2);
    res.pva.intIt = zeros(intNt,intNt,2,2);

    qp = cdat;
    qpv = permute(cat(4, cdat, cddat),[1 4 2 3]);
    qpva = permute(cat(4, cdat, cddat, cd2dat),[1 4 2 3]);

    for t1=1:intNt
        ti1 = find(time==inttime(t1));
        for t2=1:intNt
            ti2 = find(time==inttime(t2));
            if t2==t1
                continue
            end
            for stmi=1:2
                for eleci=1:2
                    otherstm = 1 + mod(stmi,2);
                    
                    % P
                    Ic1 = cmi_ggg(cstm(:,stmi), qp(:,ti1,eleci), cstm(:,otherstm), true, false, false);
                    Ic2 = cmi_ggg(cstm(:,stmi), qp(:,ti2,eleci), cstm(:,otherstm), true, false, false);
                    Ic_j = cmi_ggg( [qp(:,ti1,eleci) qp(:,ti2,eleci)], cstm(:,stmi), cstm(:,otherstm), true, false, false);
                    res.p.intItc(t1,t2,stmi,eleci) = Ic_j - Ic1 - Ic2;
                    
                    I1 = info_gg(cstm(:,stmi), qp(:,ti1,eleci), true, false, false);
                    I2 = info_gg(cstm(:,stmi), qp(:,ti2,eleci), true, false, false);
                    I_j = info_gg( [qp(:,ti1,eleci) qp(:,ti2,eleci)], cstm(:,stmi), true, false, false);
                    res.p.intIt(t1,t2,stmi,eleci) = I_j - I1 - I2;
                    
                    % PV
                    Ic1 = cmi_ggg(cstm(:,stmi), qpv(:,:,ti1,eleci), cstm(:,otherstm), true, false, false);
                    Ic2 = cmi_ggg(cstm(:,stmi), qpv(:,:,ti2,eleci), cstm(:,otherstm), true, false, false);
                    Ic_j = cmi_ggg( [qpv(:,:,ti1,eleci) qpv(:,:,ti2,eleci)], cstm(:,stmi), cstm(:,otherstm), true, false, false);
                    res.pv.intItc(t1,t2,stmi,eleci) = Ic_j - Ic1 - Ic2;
                    
                    I1 = info_gg(cstm(:,stmi), qpv(:,:,ti1,eleci), true, false, false);
                    I2 = info_gg(cstm(:,stmi), qpv(:,:,ti2,eleci), true, false, false);
                    I_j = info_gg( [qpv(:,:,ti1,eleci) qpv(:,:,ti2,eleci)], cstm(:,stmi), true, false, false);
                    res.pv.intIt(t1,t2,stmi,eleci) = I_j - I1 - I2;
                    
                    % PVA
                    Ic1 = cmi_ggg(cstm(:,stmi), qpva(:,:,ti1,eleci), cstm(:,otherstm), true, false, false);
                    Ic2 = cmi_ggg(cstm(:,stmi), qpva(:,:,ti2,eleci), cstm(:,otherstm), true, false, false);
                    Ic_j = cmi_ggg( [qpva(:,:,ti1,eleci) qpva(:,:,ti2,eleci)], cstm(:,stmi), cstm(:,otherstm), true, false, false);
                    res.pva.intItc(t1,t2,stmi,eleci) = Ic_j - Ic1 - Ic2;
                    
                    I1 = info_gg(cstm(:,stmi), qpva(:,:,ti1,eleci), true, false, false);
                    I2 = info_gg(cstm(:,stmi), qpva(:,:,ti2,eleci), true, false, false);
                    I_j = info_gg( [qpva(:,:,ti1,eleci) qpva(:,:,ti2,eleci)], cstm(:,stmi), true, false, false);
                    res.pva.intIt(t1,t2,stmi,eleci) = I_j - I1 - I2;
                end
            end
        end
    end
    
    %
    % CROSS-ELECTRODE TEMPORAL INTERACTION
    %
    xl = [0 400];
    idx = (time>xl(1)) & (time<xl(2));
    intNt = sum(idx);
    inttime = time(idx);
    % conditional
    res.p.xintItc = zeros(intNt,intNt,2);
    res.pv.xintItc = zeros(intNt,intNt,2);
    res.pva.xintItc = zeros(intNt,intNt,2);
    % direct
    res.p.xintIt = zeros(intNt,intNt,2);
    res.pv.xintIt = zeros(intNt,intNt,2);
    res.pva.xintIt = zeros(intNt,intNt,2);

    qp = cdat;
    qpv = permute(cat(4, cdat, cddat),[1 4 2 3]);
    qpva = permute(cat(4, cdat, cddat, cd2dat),[1 4 2 3]);

    for t1=1:intNt
        ti1 = find(time==inttime(t1));
        for t2=1:intNt
            ti2 = find(time==inttime(t2));
            for stmi=1:2
                otherstm = 1 + mod(stmi,2);
                
                % P
                Ic1 = cmi_ggg(cstm(:,stmi), qp(:,ti1,1), cstm(:,otherstm), true, false, false);
                Ic2 = cmi_ggg(cstm(:,stmi), qp(:,ti2,2), cstm(:,otherstm), true, false, false);
                Ic_j = cmi_ggg( [qp(:,ti1,1) qp(:,ti2,2)], cstm(:,stmi), cstm(:,otherstm), true, false, false);
                res.p.xintItc(t1,t2,stmi) = Ic_j - Ic1 - Ic2;
                
                I1 = info_gg(cstm(:,stmi), qp(:,ti1,1), true, false, false);
                I2 = info_gg(cstm(:,stmi), qp(:,ti2,2), true, false, false);
                I_j = info_gg( [qp(:,ti1,1) qp(:,ti2,2)], cstm(:,stmi), true, false, false);
                res.p.xintIt(t1,t2,stmi) = I_j - I1 - I2;
                
                % PV
                Ic1 = cmi_ggg(cstm(:,stmi), qpv(:,:,ti1,1), cstm(:,otherstm), true, false, false);
                Ic2 = cmi_ggg(cstm(:,stmi), qpv(:,:,ti2,2), cstm(:,otherstm), true, false, false);
                Ic_j = cmi_ggg( [qpv(:,:,ti1,1) qpv(:,:,ti2,2)], cstm(:,stmi), cstm(:,otherstm), true, false, false);
                res.pv.xintItc(t1,t2,stmi) = Ic_j - Ic1 - Ic2;
                
                I1 = info_gg(cstm(:,stmi), qpv(:,:,ti1,1), true, false, false);
                I2 = info_gg(cstm(:,stmi), qpv(:,:,ti2,2), true, false, false);
                I_j = info_gg( [qpv(:,:,ti1,1) qpv(:,:,ti2,2)], cstm(:,stmi), true, false, false);
                res.pv.xintIt(t1,t2,stmi) = I_j - I1 - I2;
                
                % PVA
                Ic1 = cmi_ggg(cstm(:,stmi), qpva(:,:,ti1,1), cstm(:,otherstm), true, false, false);
                Ic2 = cmi_ggg(cstm(:,stmi), qpva(:,:,ti2,2), cstm(:,otherstm), true, false, false);
                Ic_j = cmi_ggg( [qpva(:,:,ti1,1) qpva(:,:,ti2,2)], cstm(:,stmi), cstm(:,otherstm), true, false, false);
                res.pva.xintItc(t1,t2,stmi) = Ic_j - Ic1 - Ic2;
                
                I1 = info_gg(cstm(:,stmi), qpva(:,:,ti1,1), true, false, false);
                I2 = info_gg(cstm(:,stmi), qpva(:,:,ti2,2), true, false, false);
                I_j = info_gg( [qpva(:,:,ti1,1) qpva(:,:,ti2,2)], cstm(:,stmi), true, false, false);
                res.pva.xintIt(t1,t2,stmi) = I_j - I1 - I2;
            end
        end
    end
    
    %
    % CROSS-ELECTRODE TEMPORAL INTERACTION
    %
    xl = [0 400];
    idx = (time>xl(1)) & (time<xl(2));
    intNt = sum(idx);
    inttime = time(idx);

    res.p.xintIelec = zeros(intNt,intNt);
    res.pv.xintIelec = zeros(intNt,intNt);
    res.pva.xintIelec = zeros(intNt,intNt);

    qp = cdat;
    qpv = permute(cat(4, cdat, cddat),[1 4 2 3]);
    qpva = permute(cat(4, cdat, cddat, cd2dat),[1 4 2 3]);

    for t1=1:intNt
        ti1 = find(time==inttime(t1));
        for t2=1:intNt
            ti2 = find(time==inttime(t2));
            
            % P
            res.p.xintIelec(t1,t2) = info_gg(qp(:,ti1,1), qp(:,ti2,2), true, false, false);
            
            % PV
            res.pv.xintIelec(t1,t2) = info_gg(qpv(:,:,ti1,1), qpv(:,:,ti2,2), true, false, false);

            % PVA
            res.pva.xintIelec(t1,t2) = info_gg(qpva(:,:,ti1,1), qpva(:,:,ti2,2), true, false, false);
        end
    end

    
%     fname = sprintf('temporal_info_%s_%s.mat',subid,fqnames{flti});
    fname = sprintf('temporal_info_MIelec_%s_%s.mat',subid,fqnames{flti});
    savefaststruct(fullfile(data_dir,fname),res);
end

%%

%
% Run temporal interaction - direct information between time points
%

parfor subi=1:Ns
    subid = subjects{subi};
    fname = sprintf('data_%s_%s.mat',subid,fqnames{flti});
    dat = matfile(fullfile(data_dir,fname));
    eyestm = dat.eyebubs;
    otl = dat.ferp(:,:,dat.LEmaxMI);
    otr = dat.ferp(:,:,dat.REmaxMI);
    fdat = cat(3, otl, otr);
    
    ddat = zeros(size(fdat));
    d2dat = zeros(size(fdat));
    for eli=1:2
        for trli=1:size(fdat,1)
            ddat(trli,:,eli) = gradient(squeeze(fdat(trli,:,eli)),2);
            d2dat(trli,:,eli) = gradient(ddat(trli,:,eli),2);
        end
    end
    
    cstm = copnorm(eyestm);
    cdat = copnorm(fdat);
    cddat = copnorm(ddat);
    cd2dat = copnorm(d2dat);
    
    [Ntrl, Nt, ~] = size(fdat);
    
%     fname = sprintf('temporal_info_%s_%s.mat',subid,fqnames{flti});
    fname = sprintf('temporal_info_MIelec_%s_%s.mat',subid,fqnames{flti});
    res = load(fullfile(data_dir,fname));
    
    %%%%%%
    % temporal interaction - direct information between electrodes
    % needed to properly normalise redundancy
    % WITHIN-ELECTRODE TEMPORAL INTERACTION
    %
    xl = [0 400];
    idx = (time>xl(1)) & (time<xl(2));
    intNt = sum(idx);
    inttime = time(idx);

    res.p.intIelec = zeros(intNt,intNt,Nele);
    res.pv.intIelec = zeros(intNt,intNt,Nele);
    res.pva.intIelec = zeros(intNt,intNt,Nele);

    qp = cdat;
    qpv = permute(cat(4, cdat, cddat),[1 4 2 3]);
    qpva = permute(cat(4, cdat, cddat, cd2dat),[1 4 2 3]);

    for eli=1:Nele
        for t1=1:intNt
            ti1 = find(time==inttime(t1));
            for t2=1:intNt
                if t2==t1
                    continue
                end
                ti2 = find(time==inttime(t2));
                
                % P
                res.p.intIelec(t1,t2) = info_gg(qp(:,ti1,eli), qp(:,ti2,eli), true, false, false);
                
                % PV
                res.pv.intIelec(t1,t2) = info_gg(qpv(:,:,ti1,eli), qpv(:,:,ti2,eli), true, false, false);
                
                % PVA
                res.pva.intIelec(t1,t2) = info_gg(qpva(:,:,ti1,eli), qpva(:,:,ti2,eli), true, false, false);
            end
        end
    end

    savefaststruct(fullfile(data_dir,fname),res);
end

%%
% % 
% % %
% % % Run additional info calculation
% % %
% % 
% % parfor subi=1:Ns
% %     subid = subjects{subi};
% %     fname = sprintf('data_%s_%s.mat',subid,fqnames{flti});
% %     dat = matfile(fullfile(data_dir,fname));
% %     eyestm = dat.eyebubs;
% %     otl = dat.ferp(:,:,dat.LE);
% %     otr = dat.ferp(:,:,dat.RE);
% %     fdat = cat(3, otl, otr);
% %     
% %     ddat = zeros(size(fdat));
% %     d2dat = zeros(size(fdat));
% %     for eli=1:2
% %         for trli=1:size(fdat,1)
% %             ddat(trli,:,eli) = gradient(squeeze(fdat(trli,:,eli)),2);
% %             d2dat(trli,:,eli) = gradient(ddat(trli,:,eli),2);
% %         end
% %     end
% %     
% %     cstm = copnorm(eyestm);
% %     cdat = copnorm(fdat);
% %     cddat = copnorm(ddat);
% %     cd2dat = copnorm(d2dat);
% %     
% %     [Ntrl, Nt, ~] = size(fdat);
% %     
% %     fname = sprintf('temporal_info_%s_%s.mat',subid,fqnames{flti});
% %     res = load(fullfile(data_dir,fname));
% %     
% %     %%%%%%
% %     % ADD CALCULATION HERE
% % 
% %     savefaststruct(fullfile(data_dir,fname),res);
% % end