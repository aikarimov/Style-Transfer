close all
nbins = 100;
%file name
%FName = {'Levitan_flower','Seurat_flower','VanGogh_flower','Levitan_evening','Seurat_2','VanGogh_iris','VanGogh_house','VGH2'};
FName = {'VanGogh_Chestnut','Seurat_Bridge','Levitan_Oak','F752','Seurat_Bridge_grad','Levitan_Oak_grad','VanGogh_Wheat'};
%display name - add dates, correct names...
%Name = {'Levitan flower','Seurat flower','VanGogh flower','Levitan evening','Seurat Port','VanGogh iris','VanGogh house','VGH2'};
Name = {'VanGogh Chestnut','Seurat Bridge','Levitan Oak','VanGogh Chestnut grad','Seurat Bridge grad','Levitan Oak grad','VanGogh Wheat'};
nfiles = length(Name);
nx = 3; ny = 3; %supltols
fsize = 10; %font size
fsizelegend = 6; %font size for legend
lwidth = 1; %line width for auxiliary lines

for n=1:nfiles
    %xl.LabelVerticalAlignment = 'up';
    %xl.LabelHorizontalAlignment = 'left';
    load([FName{1,n},'.mat']);
    M = length(strokes);
    broadness = zeros(1,M);
    elongatedness = zeros(1,M);
    strokelengths = zeros(1,M);
    BH = zeros(1,M); %broadness homogenity
    
    for i=1:M
        broadness(1,i) = mean(strokes{1,i}.Ws);
        BH(1,i) = std(strokes{1,i}.Ws);
        elongatedness(1,i) = strokes{1,i}.length/broadness(1,i);
        strokelengths(1,i) = strokes{1,i}.length;
    end

    
    %first, get centers of brushstrokes
    strcent = zeros(2,M);
    for i=1:M
        strcent(1,i) = mean(strokes{1,i}.Xs);
        strcent(2,i) = mean(strokes{1,i}.Ys);
    end

    %then, get orientation of brushstrokes - alternative to definition in
    %Li, 2011, we obtain it from least square linear fit
    angles = zeros(1,M);
    for i = 1:M
        p = polyfit(strokes{1,i}.Xs,strokes{1,i}.Ys,1); %get linear fit
        angles(i) = atan(p(1)); %derivative = tan(angle)
    end

    %then, get straightness of brushstrokes
    straightness = zeros(1,M);
    for i = 1:M
        Np = length(strokes{1,i}.Xs);
        Suv = Np * sum(strokes{1,i}.Xs.*strokes{1,i}.Ys) - sum(strokes{1,i}.Xs)*sum(strokes{1,i}.Ys);
        Su = sqrt(Np * sum(strokes{1,i}.Xs.^2) - (sum(strokes{1,i}.Xs))^2);
        Sv = sqrt(Np * sum(strokes{1,i}.Ys.^2) - (sum(strokes{1,i}.Ys))^2);
        if (Su*Sv) == 0
            if Suv == 0
                straightness(i) = 1;
            else
                straightness(i) = abs(Suv)/1e-7;
            end
        else
            straightness(i) = abs(Suv)/(Su*Sv);
        end
    end
    
    %calculate Number of brushstrokes in the neighborhood (NBS-NB)
    %calculate Number of brushstrokes in the neighborhood with similar orientation (NBS-SO)
    %calculate orientation standard deviation in the neighborhood (OSD-NB)
    s = 200; %threshold for neighrours
    al = 0.35; %threshold for orientation
    nbsnb = zeros(1,M);
    nbsso = zeros(1,M);
    osdnb = zeros(1,M);
    for i=1:M
        nbsan = zeros(1,M); %angles of all neighbours
        for j = 1:M 
            if abs(strcent(1,i) - strcent(1,j)) < s && abs(strcent(2,i) - strcent(2,j)) < s && i~= j %criterion of being in neighbourhood
                nbsnb(i) = nbsnb(i) + 1;
                if abs(angles(i) - angles(j)) < al %criterion of similar orientation
                    nbsso(i) = nbsso(i) + 1;
                end
                nbsan(nbsnb(i)) = angles(j);
            end
        end
        osdnb(i) = std(nbsan(1:nbsnb(i))); %standard deviation of angles
    end
    
    figure(11); %нарисуем гистограмму длин
    subplot(nx,ny,n);
    histogram(strokelengths,nbins,"Normalization","probability","EdgeColor","none");
    title([Name{1,n}],'FontSize',fsize);
    xlabel('length, mm','FontSize',fsize);
    ylabel('count','FontSize',fsize);
    xlim([0 15]);
    ylim([0 0.2]);
    grid on
    

%     xline(mean(strokelengths),'-',{'Avarage'},'FontSize',fsize,'LineWidth',lwidth);
%     xline(median(strokelengths),'-',{'Median'},'FontSize',fsize,'LineWidth',lwidth);

    xline(mean(strokelengths),'--k','LineWidth',lwidth);
    xline(median(strokelengths),'-r','LineWidth',lwidth);

    legend('Data',['Average =',num2str(mean(strokelengths))], ['Median =',num2str(median(strokelengths))],'FontSize',fsizelegend)
    
    figure(12); %нарисуем гистограмму
    subplot(nx,ny,n);
    histogram(broadness,nbins,"Normalization","probability","EdgeColor","none");
    title([Name{1,n}],'FontSize',fsize);
    xlabel('width, mm','FontSize',fsize);
    ylabel('count','FontSize',fsize);
    xlim([3.9 7]);
    ylim([0 0.2]);
    grid on

    
%     xline(mean(broadness),'-',{'Avarage'},'FontSize',fsize,'LineWidth',lwidth);
%     xline(median(broadness),'-',{'Median'},'FontSize',fsize,'LineWidth',lwidth);

    xline(mean(broadness),'--k','LineWidth',lwidth);
    xline(median(broadness),'-r','LineWidth',lwidth);
    legend('Data',['Average =',num2str(mean(broadness))], ['Median =',num2str(median(broadness))],'FontSize',fsizelegend)
    
    figure(13); %нарисуем гистограмму
    subplot(nx,ny,n);
    histogram(elongatedness,nbins,"Normalization","probability","EdgeColor","none");
    title([Name{1,n}],'FontSize',fsize);
    xlabel('elongatedness, mm','FontSize',fsize);
    ylabel('count','FontSize',fsize);
    xlim([0 4]);
    ylim([0 0.16]);
    grid on
%     xline(mean(elongatedness),'-',{'Avarage'},'FontSize',fsize,'LineWidth',lwidth);
%     xline(median(elongatedness),'-',{'Median'},'FontSize',fsize,'LineWidth',lwidth);

    xline(mean(elongatedness),'--k','LineWidth',lwidth);
    xline(median(elongatedness),'-r','LineWidth',lwidth);
    legend('Data',['Average =',num2str(mean(elongatedness))], ['Median =',num2str(median(elongatedness))],'FontSize',fsizelegend)


    %Orientation
    figure(14);
    subplot(nx,ny,n);
    h = histogram(angles,nbins,"Normalization","probability","EdgeColor","none");
    polarplot(h.BinEdges([1:end-2,end]),h.Values,'LineWidth',lwidth);
    title([Name{1,n}],'FontSize',fsize);
    thetalim([-90; 90]);
    rticks([0]);
    thetaticks([-90,-45,0,45, 90]);
    thetaticklabels({'-90';'-45';'0'; '45'; '90'});
    grid on
   

    %Orientation histogram
    figure(104);
    subplot(nx,ny,n);
    histogram(angles,nbins,"Normalization","probability","EdgeColor","none");
    title([Name{1,n}],'FontSize',fsize);
    xlabel('orientation, rad','FontSize',fsize);
    ylabel('count','FontSize',fsize);
    xlim([-pi/2 pi/2]);
    ylim([0 0.1]);
    grid on
    xline(mean(angles),'--k','LineWidth',lwidth);
    xline(median(angles),'-r','LineWidth',lwidth);
    legend('Data',['Average =',num2str(mean(angles))], ['Median =',num2str(median(angles))],'FontSize',fsizelegend)


    %straightness
    figure(15);
    subplot(nx,ny,n);
    histogram(straightness,nbins,"Normalization","probability","EdgeColor","none");
    title([Name{1,n}],'FontSize',fsize);
    xlabel('straightness','FontSize',fsize);
    ylabel('count','FontSize',fsize);
    xlim([0.4 1]);
    ylim([0 0.2]);
    grid on
    xline(mean(straightness),'--k','LineWidth',lwidth);
    xline(median(straightness),'-r','LineWidth',lwidth);
    legend('Data',['Average =',num2str(mean(straightness))], ['Median =',num2str(median(straightness))],'FontSize',fsizelegend)

    %nbsnb
    figure(16);
    subplot(nx,ny,n);
    histogram(nbsnb/M,nbins,"Normalization","probability","EdgeColor","none");
    title([Name{1,n}],'FontSize',fsize);
    xlabel('NBS-NB','FontSize',fsize);
    ylabel('count','FontSize',fsize);
    xlim([0.1 0.7]);
    ylim([0 0.03]);
    grid on
    xline(mean(nbsnb),'--k','LineWidth',lwidth);
    xline(median(nbsnb),'-r','LineWidth',lwidth);
    legend('Data',['Average =',num2str(mean(nbsnb))], ['Median =',num2str(median(nbsnb))],'FontSize',fsizelegend)

    %nbsso
    figure(17);
    subplot(nx,ny,n);
    histogram(nbsso/M,nbins,"Normalization","probability","EdgeColor","none");
    title([Name{1,n}],'FontSize',fsize);
    xlabel('NBS-SO','FontSize',fsize);
    ylabel('count','FontSize',fsize);
    xlim([0 0.6]);
    ylim([0 0.03]);
    grid on
    xline(mean(nbsso),'--k','LineWidth',lwidth);
    xline(median(nbsso),'-r','LineWidth',lwidth);
    legend('Data',['Average =',num2str(mean(nbsso))], ['Median =',num2str(median(nbsso))],'FontSize',fsizelegend)

    %osdnb
    figure(18);
    subplot(nx,ny,n);
    histogram(osdnb,nbins,"Normalization","probability","EdgeColor","none");
    title([Name{1,n}],'FontSize',fsize);
    xlabel('OSD-NB','FontSize',fsize);
    ylabel('count','FontSize',fsize);
    xlim([0 1.4]);
    ylim([0 0.04]);
    grid on
    xline(mean(osdnb),'--k','LineWidth',lwidth);
    xline(median(osdnb),'-r','LineWidth',lwidth);
    legend('Data',['Average =',num2str(mean(osdnb))], ['Median =',num2str(median(osdnb))],'FontSize',fsizelegend)

    %BH
    figure(19);
    subplot(nx,ny,n);
    histogram(BH,nbins,"Normalization","probability","EdgeColor","none");
    title([Name{1,n}],'FontSize',fsize);
    xlabel('BH','FontSize',fsize);
    ylabel('count','FontSize',fsize);
%     xlim([0 1.2]);
%     ylim([0 400]);
    grid on
    xline(mean(BH),'--k','LineWidth',lwidth);
    xline(median(BH),'-r','LineWidth',lwidth);
    legend('Data',['Average =',num2str(mean(BH))], ['Median =',num2str(median(BH))],'FontSize',fsizelegend)

    %Elongatedness - Straightness plot
    figure(20);
    hold on
    scatter(mean(elongatedness),mean(straightness),'DisplayName',Name{1,n});
    title('Elongatedness vs Straightness','FontSize',fsize);
    xlabel('elongatedness','FontSize',fsize);
    ylabel('straightness','FontSize',fsize);
%     xlim([0 1.2]);
%     ylim([0 400]);
    grid on
    %xline(mean(BH),'--k','LineWidth',lwidth);
    %xline(median(BH),'-r','LineWidth',lwidth);
    legend;

end