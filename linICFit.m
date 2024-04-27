function [ic,n,resistance,offset,rsquare] = linICFit(critVolt,lowerCurr,upperCurr,lowerVolt,upperVolt,currents,voltages,bool2Plot)

%Linear fit to determine resistance & offset
[RandOffFit,S] = polyfit(currents(currents>=lowerCurr&currents<=upperCurr),voltages(currents>=lowerCurr&currents<=upperCurr),1);
[Y,DELTA] = polyconf(RandOffFit,currents,S);
resistance = RandOffFit(1);
offset = RandOffFit(2);
%Adjust data for calculated resistance and offset
rvoltages = voltages;
for i = 1:size(currents,1)
    rvoltages(i) = rvoltages(i) - ((resistance*currents(i)) + offset);
end

rcurrents = currents(rvoltages>lowerVolt&rvoltages<upperVolt);
rvoltages = rvoltages(rvoltages>lowerVolt&rvoltages<upperVolt);
while(std(rcurrents) > 10)
rcurrents = rcurrents(2:length(rcurrents));
rvoltages = rvoltages(2:length(rvoltages));
end
%Actual Ic Fitting
[icNFit,S] = polyfit(log(rcurrents),log(rvoltages),1);
[Y2,DELTA2] = polyconf(icNFit,log(currents),S);
ic = exp((log(critVolt)-icNFit(2))/icNFit(1));
n = icNFit(1);

%Calculates R-Squared for fit
resSumSquares = 0;
totalSumSquares = 0;
meanVolt = mean(voltages);
for i = 1:size(currents,1)
    resSumSquares = resSumSquares + ((voltages(i) - ((critVolt*((currents(i)/ic).^n)) + (resistance*currents(i))+offset)).^2);
    totalSumSquares = totalSumSquares + ((voltages(i) - meanVolt).^2);
end
rsquare = 1-(resSumSquares/totalSumSquares);

if(bool2Plot)
    % Plot fit with data.
    hold on;
    x = min(currents):0.1:max(currents);
    y = (critVolt*((x/ic).^n)) + (resistance.*x) + offset;
    scatter(currents,voltages,'k','filled','DisplayName','Raw Data');
    plot(x,y,'Color','r','LineWidth',2);
    xlabel('Current [A]','FontSize',20);
    ylabel('Voltage [V]','FontSize',20);
    title('Resultant Fit and Raw Data','FontSize',20);
    temp = legend;
    temp = string(temp.String);
    xhold = xlim;
    yhold = ylim;
    plot([lowerCurr,lowerCurr],[yhold(1),yhold(2)],'c','LineWidth',3,'HandleVisibility','off');
    plot([upperCurr,upperCurr],[yhold(1),yhold(2)],'c','LineWidth',3);
    plot([xhold(1),xhold(2)],[lowerVolt,lowerVolt]+offset,'m','LineWidth',3,'HandleVisibility','off');
    plot([xhold(1),xhold(2)],[upperVolt,upperVolt]+offset,'m','LineWidth',3);
    plot(currents,Y-DELTA+exp(Y2-DELTA2),'g','LineWidth',2,'HandleVisibility','off');
    plot(currents,Y+DELTA+exp(Y2+DELTA2),'g','LineWidth',2);
    plot([xhold(1),xhold(2)],[critVolt,critVolt]+offset,'b','LineWidth',2);
    legend([temp(1:length(temp)-1),sprintf(strcat('I_c: \t',num2str(ic),'\nn: \t',num2str(n),'\nResistance: \t',num2str(resistance),'\nOffset: \t',num2str(offset),'\nR^2: \t',num2str(rsquare))),'Flat Region','Transition Region','95% Confidence Interval','Critical Voltage'],'Location','northwest','FontSize',18);
    ylim([min(voltages),max(voltages)]);
    hold off;
end

end

