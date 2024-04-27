function icvb = icBFit(fields,ics,bool2plot)
    fields = round(fields,1);   %Round to 1 decimal place for fields
    uniqueFields = unique(fields);  %Find unique rounded field values
    avgIcs = zeros(length(uniqueFields),1); %Initialize array of avg Ics of length of unique fields
    for i = 1:length(uniqueFields)
        avgIcs(i) = mean(ics(uniqueFields(i) == fields));   %Average Ics at each unique field
    end
    temp = polyfit(log(abs(uniqueFields)),log(abs(avgIcs)),1);  %Linfit of log-log of Ic v B
    icvb.alpha = temp(1);   %Fill output variable
    icvb.ic0 = exp(temp(2));    %Fill output variable
    if(bool2plot)
            hold on;
            xlabel('Magnetic Field [T]','FontSize',25);
            ylabel('I_c [A]','FontSize',25);
            title('I_c vs. B');
            scatter(fields,ics,'filled');
            plot(min(uniqueFields):0.1:max(uniqueFields),icvb.ic0*((min(uniqueFields):0.1:max(uniqueFields)).^icvb.alpha),'HandleVisibility','off','LineWidth',1.75);
            temp = legend;
            temp = temp.String;
            legend([temp(1:length(temp)-1),sprintf(strcat('Alpha: \t',num2str(-1*icvb.alpha),'\nI_c_0: \t',num2str(icvb.ic0),'\t'))],'FontSize',20);
            hold off;
            %Dynamically adds to prexisting graph
     end
end

