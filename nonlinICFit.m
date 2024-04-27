function [ic,n,resistance,offset,rsquare] = nonlinICFit(critVolt,currents,voltages,bool2Plot)
if (length(voltages) > 10)
    % Define custom non-linear regression equation to approach
    ft = fittype('((x/a)^b) + c*x + d','independent','x','dependent','y');
    opts = fitoptions('Method','NonlinearLeastSquares');
    opts.DiffMinChange = 1e-4;
    opts.Display = 'Off';
    opts.MaxFunEvals = 1000;
    opts.MaxIter = 1000;
    opts.TolFun = 1e-6;
    opts.TolX = 1e-6;
    [xData,yData] = prepareCurveData(currents,voltages/critVolt);
    yData = yData(xData>5);
    xData = xData(xData>5);
    %Remember, y-values being fitting are actually V/Vc
    opts.Lower = [1e-6 0 -(1e-4)/critVolt -(1e-2)/critVolt];
    
    %"Reasonable" initial guesses for coefficients
    flipCurr = currents';
    flipVolt = voltages';
    holdMe = 1:1:floor(length(flipCurr)/2);
    Roffguess = polyfit(flipCurr(holdMe),flipVolt(holdMe),1);
    opts.StartPoint = [max(currents)*0.9 15 Roffguess(1)/critVolt Roffguess(2)/critVolt];
    
    %Remember, y-values being fitted are actually V/Vc
    opts.Upper = [max(currents) 200 (1e-4)/critVolt (1e-2)/critVolt];
    %All values above were guesses unless otherwise stated
    
    
    % Fit model to data.
    try
        [fitresult, gof] = fit(xData,yData,ft,opts);
        bounds = confint(fitresult);
        
        %Put everything interesting into outputs
        coeffs = coeffvalues(fitresult);
        ic = coeffs(1);
        n = coeffs(2);
        resistance = coeffs(3)*critVolt;
        offset = coeffs(4) * critVolt;
        rsquare = gof.rsquare;
        bounds(1,3) = bounds(1,3) * critVolt;
        bounds(1,4) = bounds(1,4) * critVolt;
        bounds(2,3) = bounds(2,3) * critVolt;
        bounds(2,4) = bounds(2,4) * critVolt;
    catch
        ic = 0;
        n = 0;
        resistance = 0;
        offset = 0;
        rsquare = 0;
    end
    if(bool2Plot)
        % Plot fit with data.
        figure('Name','IV Curve Fit');
        hold on;
        scatter(currents,voltages,[],'k','*','DisplayName','Raw Data');
        I = min(currents):0.01:max(currents);
        I = I(I>0);
        plot(I,(critVolt*((I/ic).^n))+(resistance*I)+offset,'r','LineWidth',2,'DisplayName',sprintf(strcat('I_c: \t',num2str(round(ic,1)),'\nn: \t',num2str(round(n,1)),'\nResistance: \t',num2str(round(resistance,11)),'\nOffset: \t',num2str(round(offset,7)),'\nR^2: \t',num2str(round(rsquare,3)))));
        xlabel('Current [A]','FontSize',35);
        ylabel('Voltage [V]','FontSize',35);
        title('Resultant Fit and Raw Data','FontSize',40);
        plot(I,(critVolt*((I/bounds(1,1)).^bounds(2,2)))+(bounds(2,3)*I)+bounds(2,4),'g','LineWidth',2,'HandleVisibility','off');
        plot(I,(critVolt*((I/bounds(2,1)).^bounds(1,2)))+(bounds(1,3)*I)+bounds(1,4),'g','LineWidth',2,'DisplayName','95% Confidence Interval');
        xhold = xlim;
        plot([xhold(1),xhold(2)],[critVolt,critVolt]+offset+(resistance*ic),'b','LineWidth',2,'DisplayName','Critical Voltage');
        legend('Location','northwest','FontSize',30);
        ylim([min(voltages),max(voltages)]);
        movegui('center');
        grid on;
        hold off;
    end
else
    ic = 0;
    n = 0;
    resistance = 0;
    offset = 0;
    rsquare = 0;
end

end