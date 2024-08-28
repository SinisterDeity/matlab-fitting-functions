function coefficients = icBThetaFit(obj,bool2Plot)
%it works bro

coefficients(length(obj)) = struct('a',0,'b',0,'c',0,'d',0,'e',0,'f',0,'g',0,'h',0,'i',0);
for z = length(obj)
    angles = obj(z).angles;
    fields = obj(z).fields;
    ics = obj(z).ics;

    [xData, yData, zData] = prepareSurfaceData( angles, fields, ics );
    ft = fittype( 'c*((y^((l-(m*k)/((cosd(x)^2)+((k^2)*sind(x)^2)))))*(((a*x+q)*b/((cosd(x)^2)+((b^2)*(sind(x)^2)))))+(h*p/((cosd(x)^2)+((p^2)*(sind(x)^2)))))', 'independent', {'x', 'y'}, 'dependent', 'z' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.DiffMaxChange = 10;
    opts.DiffMinChange = 1e-12;
    opts.Display = 'Off';
    opts.Lower = [-1 0.001 0.001 0.001 0 0.001 0.001 0 0];
    opts.MaxFunEvals = 10000;
    opts.MaxIter = 5000;
    opts.Robust = 'Bisquare';
    opts.StartPoint = [-1 1 1 1 1 1 1 1 1];
    opts.TolFun = 1e-12;
    opts.TolX = 1e-09;
    opts.Upper = [0 Inf Inf Inf Inf Inf Inf Inf Inf];


    [fitresult, gof] = fit( [xData, yData], zData, ft, opts );
    coeffs = coeffvalues(fitresult);
    coefficients(z).a = coeffs(3);
    coefficients(z).b = coeffs(6);
    coefficients(z).c = coeffs(7);
    coefficients(z).d = coeffs(5);
    coefficients(z).e = coeffs(1);
    coefficients(z).f = coeffs(9);
    coefficients(z).g = coeffs(2);
    coefficients(z).h = coeffs(4);
    coefficients(z).i = coeffs(8);

    if(bool2Plot)
        hold on;
        scatter3(angles,fields,ics,'k','filled');
        surf(angles,fields,icBThetaEval(coefficients,angles,fields));
        xlabel('Angle [°]');
        ylabel('Magnetic Field [T]');
        zlabel('I_c [A]');
        title('$I_c(B (15 - 30T),\theta (-20,20), T (4.2 K))= a \left( \left(\left(B^{b - \frac{c \cdot d}{\cos^2(\theta) + d^2 \cdot \sin^2(\theta)}}\right) \left( \frac{(e \cdot B + f) \cdot g}{\cos^2(\theta) + g^2 \cdot \sin^2(\theta)} \right)\right) + \frac{h \cdot k}{\cos^2(\theta) + k^2 \cdot \sin^2(\theta)} \right)$','Interpreter','latex');
        view( 183.1, -0.6 );
        grid on;
        ax = gca;
        ax.FontSize=20;
    end
end