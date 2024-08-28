function ic = singleLorentz(obj,angle)
%Input struct must contain .i0,.gamma

ic = zeros(length(angle),length(obj));
for j = 1:length(obj)
    i0 = obj(j).i0;
    gamma = obj(j).gamma;
    for i =1:length(angle)
        ic(i,j) = ((i0.*gamma./pi).*(((cosd(angle(i)).^2)+((gamma.^2).*(sind(angle(i)).^2))).^-1));
    end
end
end