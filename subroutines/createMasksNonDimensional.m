function [topmat, botmat] = createMasksNonDimensional(Z,eta,epsilon,delta)



[mm nn] = size(Z);
botmat = ones(mm,nn);
topmat = botmat;
for m = 1:mm
    for n = 1:nn
        if (Z(m,n)>epsilon*eta(n))
            botmat(m,n)=NaN;
        elseif Z(m,n)<epsilon*eta(n)
            topmat(m,n) = NaN;
            
        end
    end
end
