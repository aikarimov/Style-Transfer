close all;

hval = [0 0 0 0 0.31 0.31 1.23 1.23 0.32 0.65 0.65 1.23 1.86 1.86 1.86 2.41 3.67 5.77 7.29 10.81 11.54 13.54 11.54 17.79 15.79 17.79 20.96 13.61 17.6 14.82 15.96 13.59 16.64 15.45 23.89 28.48 38.77 42.72 38.78 32.98 22.09 13.26 6.28 1.16 0];
n = length(hval);
edges = linspace(0,180,n+1);
cf = 0.082/max(hval);
hval = hval*cf; %scaling
%make syntetic dataset
nstrokes = 2500; %estimated
nbsnbr = []; %real nbsnb
for i = 1:n
    nbsnbr = [nbsnbr,repmat(0.5*(edges(i) + edges(i+1)),1,round(nstrokes*hval(i)))];
end
%normalization
nbsnbr = nbsnbr / nstrokes;

nbins = n;
%file name
FName = {'VanGogh_Chestnut','VanGogh_Chestnut_grad'};
%display name - add dates, correct names...
Name = {'Van Gogh (NST)','Van Gogh (GRAD)','Van Gogh (Real)'};
nfiles = length(FName);
nx = 1; ny = 3; %supltols
fsize = 10; %font size
fsizelegend = 6; %font size for legend
lwidth = 1; %line width for auxiliary lines

%create matrices for correlation analysis
cnbsnb = zeros(nbins,nfiles);

for n=1:nfiles
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
   
    %calculate Number of brushstrokes in the neighborhood (NBS-NB)
    s = 200; %threshold for neighrours
    al = 0.35; %threshold for orientation
    nbsnb = zeros(1,M);
    for i=1:M
        for j = 1:M 
            if abs(strcent(1,i) - strcent(1,j)) < s && abs(strcent(2,i) - strcent(2,j)) < s && i~= j %criterion of being in neighbourhood
                nbsnb(i) = nbsnb(i) + 1;
            end
        end
    end

    %normalization
    nbsnb = nbsnb/M;

    %nbsnb
    figure(16);
    subplot(nx,ny,n);
    h = histogram(nbsnb,nbins,"Normalization","probability","EdgeColor","none");
    cnbsnb(:,n) = (h.Values)';
    title([Name{1,n}],'FontSize',fsize);
    xlabel('NBS-NB','FontSize',fsize);
    ylabel('count','FontSize',fsize);
    xlim([0.1 0.7]);
    ylim([0 0.1]);
    grid on
    xline(mean(nbsnb),'--k','LineWidth',lwidth);
    xline(median(nbsnb),'-r','LineWidth',lwidth);
    legend('Data',['Average =',num2str(mean(nbsnb))], ['Median =',num2str(median(nbsnb))],'FontSize',fsizelegend)
end

%plot real
%nbsnb
n = 3;
figure(16);
subplot(nx,ny,n);
h = histogram('BinEdges',edges/327,'BinCounts',hval,"EdgeColor","none");
cnbsnb(:,n) = (h.Values)';
title([Name{1,n}],'FontSize',fsize);
xlabel('NBS-NB','FontSize',fsize);
ylabel('count','FontSize',fsize);
xlim([0.1 0.7]);
ylim([0 0.1]);
grid on
xline(mean(nbsnbr),'--k','LineWidth',lwidth);
xline(median(nbsnbr),'-r','LineWidth',lwidth);
legend('Data',['Average =',num2str(mean(nbsnbr))], ['Median =',num2str(median(nbsnbr))],'FontSize',fsizelegend)


figure(24);
cdata = corrcoef(cnbsnb);
h = heatmap(Name,Name,cdata);
h.Title = 'NBS-NB histogram correlation';

