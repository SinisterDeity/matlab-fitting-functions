function result = icBThetaEval(obj,angles,fields)
%obj must feed in a,b,c,d,e,f,g,h,i as struct params

a=obj.a;
b=obj.b;
c=obj.c;
d=obj.d;
e=obj.e;
f=obj.f;
g=obj.g;
h=obj.h;
i=obj.i;
result=a.*((fields.^((b-(c.*d)./((cosd(angles).^2)+((d.^2)*sind(angles).^2))))).*(((e.*angles+f).*g./((cosd(angles).^2)+((g.^2).*(sind(angles).^2)))))+(h.*i./((cosd(angles).^2)+((i.^2)*(sind(angles).^2)))));
end