function [topmat, botmat] = createMasks(Z,eta,h1,h2)



[mm nn] = size(Z);
botmat = ones(mm,nn);
topmat = botmat;
for m = 1:mm
    for n = 1:nn
        if (Z(m,n)>eta(n))||(Z(m,n)<-h2)
            botmat(m,n)=NaN;
        end
        if (Z(m,n)<eta(n))||(Z(m,n)>h1)
            topmat(m,n) = NaN;
        end
    end
end
