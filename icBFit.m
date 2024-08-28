function ic0sAlphas = icBFit(obj,bool2plot)
%ICBFIT

parfor z = 1:length(obj)
    holder = obj(z);
    fields = holder.fields';
    ics = holder.ics';
    fields = round(fields,1);   %Round to 1 decimal place for fields
    uniqueFields = unique(fields);  %Find unique rounded field values
    avgIcs = zeros(length(uniqueFields),1); %Initialize array of avg Ics of length of unique fields
    for i = 1:length(uniqueFields)
        avgIcs(i) = mean(ics(uniqueFields(i) == fields));   %Average Ics at each unique field
    end
    [xData, yData] = prepareCurveData(uniqueFields,avgIcs);
    ft = fittype( 'a*x^-b', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.DiffMaxChange = 1;
    opts.DiffMinChange = 1e-09;
    opts.Display = 'Off';
    opts.Lower = [0 0];
    opts.MaxFunEvals = 1000;
    opts.MaxIter = 500;
    opts.Robust = 'Bisquare';
    opts.StartPoint = [1000 0.5];
    opts.TolFun = 1e-09;
    opts.TolX = 1e-09;
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    bounds = confint(fitresult);
    coeffs = coeffvalues(fitresult);
    icvb(z).ic0 = coeffs(1);
    icvb(z).alpha = coeffs(2);
   if(bool2plot)
    hold on;
    title('I_c(B) = I_{c_{0}}(B^{-\alpha})','FontSize',25);
    xlabel('Magnetic Field [T]','FontSize',25);
    ylabel('I_c [A]','FontSize',25);
    scatter(obj(z).fields,obj(z).ics,'k','filled','DisplayName','Raw Data');
    fields = min(obj(z).fields):0.1:max(obj(z).fields);
    ics = coeffs(1).*(fields.^-coeffs(2));
    plot(fields,ics,'r','LineWidth',2,'DisplayName',sprintf(strcat('I_{c_{0}} = ',num2str(round(coeffs(1),1)),'\nAlpha = ',num2str(round(coeffs(2),3)),'\nR-Squared = ',num2str(round(gof.rsquare,3)))));
    ics = bounds(1,1).*((fields).^-bounds(2,2));
    plot(fields,ics,'g','LineWidth',2,'DisplayName','Lower 95% CI');
    ics = bounds(2,1).*((fields).^-bounds(1,2));
    plot(fields,ics,'b','LineWidth',2,'DisplayName','Upper 95% CI');
    ax = gca;
    ax.FontSize = 25;
    legend;
   end
end
ic0sAlphas = [[icvb.ic0];[icvb.alpha]]';
end

