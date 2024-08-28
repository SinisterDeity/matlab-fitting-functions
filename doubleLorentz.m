function ic = doubleLorentz(obj,angle)
%Input struct must contain .i01,.gamma1,.i02,.gamma2

ic = zeros(length(angle),length(obj));
for j = 1:length(obj)
    I01 = obj(j).i01;
    Gamma1 = obj(j).gamma1;
    I02 = obj(j).i02;
    Gamma2 = obj(j).gamma2;
    for i =1:length(angle)
        ic(i,j) = ((I01.*Gamma1./pi).*(((cosd(angle(i)).^2)+((Gamma1.^2).*(sind(angle(i)).^2))).^-1))+((I02.*Gamma2./pi).*(((cosd(angle(i)).^2)+((Gamma2.^2).*(sind(angle(i)).^2))).^-1));
    end
end
end
