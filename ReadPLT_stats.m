%������ ��� �����, ������� ������. ��������� ���� ����
warning off

pixScle = 1; %���������� �������� �� -> �������
antialiasing = 0; %change to 1 if needed

[filename,pathname] = uigetfile('*.plt','�������� ����');
fid = fopen([pathname,filename],'r');
%filename = 'levitan_flower.plt';

fid = fopen(filename,'r');

c = fscanf(fid,'%s');
cmds = split(c,';'); %��������� � �� ������ �MDS: CMDS[i] - i-� �������
color = zeros(1,1,3); %��������� ����
mixtype = 1;
ptX = 0; ptY = 0;

brushWidth_mm = 2; sideXmm = 515; sideYmm = 675;

scalePLT = 1/40;%�������������� ������ PLT

HeigthY = sideYmm*40*scalePLT;

brushSize = brushWidth_mm * pixScle;%������� �����, � ��������
bs2 = ceil(brushSize/2);
bsQuad = bs2^2;
bsQuad2 = 2;

bsQuad_1 = (bs2 - 1)^2;
bsQuadp1 = (bs2 + 1)^2;

%PWscl = 300/25.4; %11.8, standard
PWscl = 1;

hw = waitbar(0,'���������� ������� PLT'); %���������� ����������, �������
N = length(cmds);
strokectr = 0;
colctr = 0; %color counter
codecolor = [0 0 0];

strokelengths = [];
curlength = 0;

strokes = {}; %empty cell array
stroke = [];

for j=1:N %�� ������ ������
    s = cmds{j,1}; %s - ��������� �������
    if(~isempty(s)) %���� ��� �� ������
        cmd = s(1:2); %����� ������ 2 �������
        switch cmd %� ����������� �� ��� ������ ��������
            case 'PU'    %������� �����, �������� � ���������� - �������� ptX, ptY ��� ����� ����������            
                i = strfind(s,',');
                sX = s(3: i - 1); 
                sY = s(i:end);
                %���������� ����� � �����
                ptX = str2double(sY);
                ptY = str2double(sX);
                %�������������� ��������
                ptX = (HeigthY*pixScle - round(ptX*scalePLT*pixScle));
                ptY = round(ptY*scalePLT*pixScle);
                
                %���� ��� ���� �����, ������� ��� � ������ strokes
                if ~isempty(stroke)
                    strokes{end + 1} = stroke;
                    strokelengths = [strokelengths; stroke.length];
                end
                stroke = struct('color',color,'Xs',ptX,'Ys',ptY,'length',0,'Ws',brushWidth_mm); %structure for the stroke

            case 'PD'
                i = strfind(s,',');
                sX = s(3: i - 1);
                sY = s(i:end); 
                pX = str2double(sY);
                pY = str2double(sX);
                
                pX = (HeigthY*pixScle - round(pX*scalePLT*pixScle));
                pY = round(pY*scalePLT*pixScle);
                
                %������� ����� ����� � ��������� �����
                stroke.length = stroke.length + sqrt((ptX - pX)^2 + (ptY - pY)^2)/pixScle;
                stroke.Xs = [stroke.Xs, pX];
                stroke.Ys = [stroke.Ys, pY];
                stroke.Ws = [stroke.Ws, brushWidth_mm];

                ptX = pX; ptY = pY; %�������� ������� ����� ��� ����������
                
           case 'PC' %����� �����
                waitbar(j/N,hw,['���������� ������� PLT ',num2str(j)]);
                i = strfind(s,',');
                spl = split(s(3:end),',');
                color(1,1,1)= str2double(spl{2,1});
                color(1,1,2) = str2double(spl{3,1});
                color(1,1,3) = str2double(spl{4,1});
                
            case 'PW' %����� ������ �����
                i = strfind(s,',');
                brushSize_txt = s(3:end);
                brushWidth_mm = str2double(brushSize_txt) * PWscl;
                brushSize = brushWidth_mm * pixScle;%������� �����
                bs2 = ceil(brushSize/2);
                bsQuad = bs2^2;
                bsQuad2 = 2;
        end
    end
end
close(hw); %������� �������
fclose(fid); %������� ����

figure; %�������� �����������
histogram(strokelengths,100);
title(['Stroke length histogram. Count = ',num2str(length(strokes))]);
xlabel('length, mm');
ylabel('count');

save([filename(1:end-4),'.mat'],'strokes');