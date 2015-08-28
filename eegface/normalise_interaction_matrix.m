function Iintn = normalise_interaction_matrix(Iint,I1,I2,II,thresh)

Nt = size(Iint,1);
Iintn = zeros(size(Iint));
Ired = -Iint;
% I = diag(Ired);
Ired(Ired<0) = 0;

for t1=1:Nt
    for t2=1:Nt
        Iintn(t1,t2) = 100*Ired(t1,t2) ./ min([I1(t1),I2(t2),II(t1,t2)]);
    end
end

% apply threshold
for ti=1:Nt
    if (I1(ti)<thresh)
        Iintn(ti,:) = 0;
    end
    if (I2(ti)<thresh)
        Iintn(:,ti) = 0;
    end
end

% Iintn = Iintn + Iintn';