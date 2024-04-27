function [tcfit] = tcFit(dankStruct,plotDesired)
for i = 1:length(dankStruct)
    %Load table data into suitable data structures
    holder = dankStruct(i);
    temperatures = holder.t;
    Moments = holder.m;
    %Normalize Moments from [0,1]
    Moments = normalize(medfilt1(Moments),'Range');

    %Filter out spikes using 3rd-order median filter
    Moments = medfilt1(Moments);

    %Calculates linear fit of first 10 points
    lowFit = polyfit(temperatures(1:10),Moments(1:10),1);

    %Calculates linear fit between 40% & 60% of normalized moment
    midPointLow = sum(Moments<0.4);
    midPointHigh = sum(Moments<0.6);
    midFit = polyfit(temperatures(midPointLow:midPointHigh),Moments(midPointLow:midPointHigh),1);

    %Calculates linear fit of last 10 points
    highPointLow = size(temperatures,1)-10;
    highPointHigh = size(temperatures,1);
    highFit = polyfit(temperatures(highPointLow:highPointHigh),Moments(highPointLow:highPointHigh),1);

    %Calculates the transition breadth by doing T@90% - T@10%
    lowerTransBound = temperatures(sum(Moments<=0.1));
    upperTransBound = temperatures(sum(Moments<=0.9));
    tcfit.breadth = upperTransBound - lowerTransBound;

    %genericTc = way people typically do it, customTc = Griffin's way
    tcfit.genericTc = (midFit(2) - highFit(2))/(highFit(1) - midFit(1));
    tcfit.customTc = (midFit(2) - lowFit(2))/(lowFit(1) - midFit(1));

    if(plotDesired)
        hold on;
        title('T_c Fitting','FontSize',25);
        xlabel('Temperature [K]','FontSize',25);
        ylabel('Normalized Moment','FontSize',25);
        plot(temperatures,Moments,'-*','DisplayName',sprintf(strcat('T_c: \t',num2str(tcfit.genericTc),'\nBreadth: \t',num2str(tcfit.breadth))),'LineWidth',2,'Color',Colors(i));
        ax = gca;
        ax.FontSize = 25;
        legend('FontSize',25);
        hold off;
    end
end
end

