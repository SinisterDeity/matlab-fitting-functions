function ic0sTstars = icTFit(obj,bool2plot)
%ICBFIT

for z = 1:length(obj)
    holder = obj(z);
    temperatures = holder.temperature;
    ics = holder.ic;
    temperatures = round(temperatures,1);   %Round to 1 decimal place for temperatures
    uniqueTemperatures = unique(temperatures);  %Find unique rounded temperature values
    avgIcs = zeros(length(uniqueTemperatures),1); %Initialize array of avg Ics of length of unique temperatures
    for i = 1:length(uniqueTemperatures)
        avgIcs(i) = mean(ics(uniqueTemperatures(i) == temperatures));   %Average Ics at each unique temperature
    end
    [xData, yData] = prepareCurveData(uniqueTemperatures',avgIcs);
    ft = fittype( 'a*exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.DiffMaxChange = 1;
    opts.DiffMinChange = 1e-09;
    opts.Display = 'Off';
    opts.Lower = [0 0];
    opts.MaxFunEvals = 1000;
    opts.MaxIter = 5000;
    opts.Robust = 'Bisquare';
    opts.StartPoint = [100 20];
    opts.TolFun = 1e-09;
    opts.TolX = 1e-09;

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    %bounds = confint(fitresult);
    coeffs = coeffvalues(fitresult);
    icvt(z).tstar = -coeffs(2);   %Fill output variable
    icvt(z).ic0 = coeffs(1);    %Fill output variable
    if(bool2plot)
        hold on;
        xlabel('Temperature [K]','FontSize',25);
        ylabel('I_c [A]','FontSize',25);
        title('I_c(T) = I_{c_{T=0}}e^{(-T/T_*)}');
        scatter(obj(z).temperature,obj(z).ic,'k','filled','DisplayName','Raw Data');
        temps = min(uniqueTemperatures):0.1:max(uniqueTemperatures);
        plot(temps,icvt(z).ic0*(exp((temps)./icvt(z).tstar)),'r','LineWidth',2,'DisplayName',sprintf(strcat('T^*: \t',num2str(round(-1*icvt(z).tstar,1)),'\nI_c_0: \t',num2str(round(icvt(z).ic0,1)),'\nR^2 = ',num2str(round(gof.rsquare,3)))));
        ics = bounds(1,1).*exp(-temps./bounds(1,2));
        plot(temps,ics,'g','LineWidth',2,'DisplayName','Lower 95% CI');
        ics = bounds(2,1).*exp(-temps./bounds(2,2));
        plot(temps,ics,'b','LineWidth',2,'DisplayName','Upper 95% CI');
        legend;
        ax = gca;
        ax.FontSize = 25;
        %Dynamically adds to prexisting graph
    end
end
ic0sTstars = [[icvt.ic0];-[icvt.tstar]]';
end

