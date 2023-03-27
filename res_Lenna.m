close all
nbins = 100;
%file name
FName = {'VanGogh_Lenna_Iris','Seurat_Lenna_Bridge','Levitan_Lenna_Oak','Lenna_grad'};
%display name - add dates, correct names...
Name = {'Lenna (Van Gogh)','Lenna (Seurat)','Lenna (Levitan)','Lenna (Gradient)'};
nfiles = length(Name);
nx = 2; ny = 2; %supltols
fsize = 10; %font size
fsizelegend = 6; %font size for legend
lwidth = 1; %line width for auxiliary lines

markers = {'^','square','o','diamond'};

%create matrices for correlation analysis

clength = zeros(nbins,nfiles);
cstraightness = zeros(nbins,nfiles);
corientation = zeros(nbins,nfiles);
cnbsnb = zeros(nbins,nfiles);
cnbsso = zeros(nbins,nfiles);
cosdnb = zeros(nbins,nfiles);

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

    %normalization
    nbsnb = nbsnb/M;
    nbsso = nbsso/M;

    
    figure(11); %нарисуем гистограмму длин
    subplot(nx,ny,n);
    h = histogram(strokelengths,nbins,"Normalization","probability","EdgeColor","none");
    clength(:,n) = (h.Values)';
    title([Name{1,n}],'FontSize',fsize);
    xlabel('length, mm','FontSize',fsize);
    ylabel('count','FontSize',fsize);
    xlim([0 15]);
    ylim([0 0.2]);
    grid on
    xline(mean(strokelengths),'--k','LineWidth',lwidth);
    xline(median(strokelengths),'-r','LineWidth',lwidth);

    legend('Data',['Average =',num2str(mean(strokelengths))], ['Median =',num2str(median(strokelengths))],'FontSize',fsizelegend)
    
%     figure(12); %нарисуем гистограмму
%     subplot(nx,ny,n);
%     histogram(broadness,nbins,"Normalization","probability","EdgeColor","none");
%     title([Name{1,n}],'FontSize',fsize);
%     xlabel('width, mm','FontSize',fsize);
%     ylabel('count','FontSize',fsize);
%     xlim([3.9 7]);
%     ylim([0 0.2]);
%     grid on
%     xline(mean(broadness),'--k','LineWidth',lwidth);
%     xline(median(broadness),'-r','LineWidth',lwidth);
%     legend('Data',['Average =',num2str(mean(broadness))], ['Median =',num2str(median(broadness))],'FontSize',fsizelegend)
    
%     figure(13); %нарисуем гистограмму
%     subplot(nx,ny,n);
%     histogram(elongatedness,nbins,"Normalization","probability","EdgeColor","none");
%     title([Name{1,n}],'FontSize',fsize);
%     xlabel('elongatedness, mm','FontSize',fsize);
%     ylabel('count','FontSize',fsize);
%     xlim([0 4]);
%     ylim([0 0.16]);
%     grid on
%     xline(mean(elongatedness),'--k','LineWidth',lwidth);
%     xline(median(elongatedness),'-r','LineWidth',lwidth);
%     legend('Data',['Average =',num2str(mean(elongatedness))], ['Median =',num2str(median(elongatedness))],'FontSize',fsizelegend)


    %Orientation - angles
%     figure(14);
%     subplot(nx,ny,n);
%     h = histogram(angles,nbins,"Normalization","probability","EdgeColor","none");
%     polarplot(h.BinEdges([1:end-2,end]),h.Values,'LineWidth',lwidth);
%     title([Name{1,n}],'FontSize',fsize);
%     thetalim([-90; 90]);
%     rticks([0]);
%     thetaticks([-90,-45,0,45, 90]);
%     thetaticklabels({'-90';'-45';'0'; '45'; '90'});
%     grid on
   

    %Orientation histogram
    figure(104);
    subplot(nx,ny,n);
    h = histogram(angles,nbins,"Normalization","probability","EdgeColor","none");
    corientation(:,n) = (h.Values)';
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
    h = histogram(straightness,nbins,"Normalization","probability","EdgeColor","none");
    cstraightness(:,n) = (h.Values)';
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
    h = histogram(nbsnb,nbins,"Normalization","probability","EdgeColor","none");
    cnbsnb(:,n) = (h.Values)';
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
    h = histogram(nbsso,nbins,"Normalization","probability","EdgeColor","none");
    cnbsso(:,n) = (h.Values)';
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
    h = histogram(osdnb,nbins,"Normalization","probability","EdgeColor","none");
    cosdnb(:,n) = (h.Values)';
    title([Name{1,n}],'FontSize',fsize);
    xlabel('OSD-NB','FontSize',fsize);
    ylabel('count','FontSize',fsize);
    xlim([0 1.4]);
    ylim([0 0.04]);
    grid on
    xline(mean(osdnb),'--k','LineWidth',lwidth);
    xline(median(osdnb),'-r','LineWidth',lwidth);
    legend('Data',['Average =',num2str(mean(osdnb))], ['Median =',num2str(median(osdnb))],'FontSize',fsizelegend)

    %BH - useless
%     figure(19);
%     subplot(nx,ny,n);
%     histogram(BH,nbins,"Normalization","probability","EdgeColor","none");
%     title([Name{1,n}],'FontSize',fsize);
%     xlabel('BH','FontSize',fsize);
%     ylabel('count','FontSize',fsize);
%     grid on
%     xline(mean(BH),'--k','LineWidth',lwidth);
%     xline(median(BH),'-r','LineWidth',lwidth);
%     legend('Data',['Average =',num2str(mean(BH))], ['Median =',num2str(median(BH))],'FontSize',fsizelegend)

    %Elongatedness - Straightness plot
    figure(20);
    hold on
    scatter(mean(nbsso),mean(straightness),'DisplayName',Name{1,n},'Marker',markers{n});
    title('NBS-SO vs Straightness','FontSize',fsize);
    xlabel('NBS-SO','FontSize',fsize);
    ylabel('straightness','FontSize',fsize);
    grid on
    legend;   
end

%plot correlation heatmaps
figure(21);
cdata = corrcoef(clength);
h = heatmap(Name,Name,cdata);
h.Title = 'length histogram correlation';

figure(22);
cdata = corrcoef(cstraightness);
h = heatmap(Name,Name,cdata);
h.Title = 'straightness histogram correlation';

figure(23);
cdata = corrcoef(corientation);
h = heatmap(Name,Name,cdata);
h.Title = 'orientation histogram correlation';

figure(24);
cdata = corrcoef(cnbsnb);
h = heatmap(Name,Name,cdata);
h.Title = 'NBS-NB histogram correlation';

figure(25);
cdata = corrcoef(cnbsso);
h = heatmap(Name,Name,cdata);
h.Title = 'NBS-SO histogram correlation';

figure(26);
cdata = corrcoef(cosdnb);
h = heatmap(Name,Name,cdata);
h.Title = 'OSD-NB histogram correlation';

