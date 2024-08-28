function offsetScaleGamma = alphaAngleFit(obj,bool2Plot)

%if(isscalar(obj))
%else
    offsetScaleGamma = zeros(length(obj),3);
    for i = 1:length(obj)
        [xData,yData] = prepareCurveData(obj(i).angle,obj(i).alpha);
        ft = fittype( '(a-(b*c)/(cosd(x)^2+(c^2)*sind(x)^2))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.DiffMaxChange = 1;
        opts.DiffMinChange = 1e-09;
        opts.Display = 'Off';
        opts.Lower = [0 0 0];
        opts.MaxFunEvals = 1000;
        opts.MaxIter = 500;
        opts.Robust = 'Bisquare';
        opts.StartPoint = [0.823457828327293 0.694828622975817 0.317099480060861];
        opts.TolFun = 1e-09;
        opts.TolX = 1e-09;

        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts );
        coeffs = coeffvalues(fitresult);
        offsetScaleGamma(i,:) = [coeffs(1),coeffs(2),coeffs(3)];

        if(bool2Plot)
            figure;
        hold on;
        title('$\alpha{}(\theta{})$','Interpreter','latex');
        xlabel('Angle [°]');
        ylabel('\alpha');
        plot(xData,invertedOffsetLorentz(offsetScaleGamma(i,:),xData));
        end
    end
%end