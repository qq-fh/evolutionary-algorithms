%% SCH
clc,clear
clear all;

data1 = load('sch1.txt');
data2 = load('sch2.txt');
data = [data1;data2];

plot(data(:, 1), data(:,2), 'ro')
axis( [-0.2, 4.5, -0.2, 4.5] )
xlabel('f_1')
ylabel('f_2')
legend('SCH')
set(gca,'FontSize',16); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(get(gca,'XLabel'),'FontSize',16);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',16);

%% FON
clc,clear
clear all;

data1 = load('fon1.txt');
data2 = load('fon2.txt');
data = [data1;data2];

% subplot(121);
% plot(data1(:, 1), data1(:,2), 'r*')
% subplot(122);
% plot(data2(:, 1), data2(:,2), 'r*')
plot(data2(:, 1), data2(:,2), 'ro')

axis( [-0.1, 1.1, -0.1, 1.1] )
xlabel('f_1')
ylabel('f_2')
legend('FON')
set(gca,'FontSize',16); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(get(gca,'XLabel'),'FontSize',16);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',16);


%% KUR
clc,clear
clear all;

data = load('kur.txt');

plot(data(:, 1), data(:,2), 'ro')

axis( [-20, -14, -12, 2] )
xlabel('f_1')
ylabel('f_2')
legend('KUR')
set(gca,'FontSize',16); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(get(gca,'XLabel'),'FontSize',16);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',16);

%% ZDT1
clc,clear
clear all;

data = load('ZDT1.txt');    %nsga2584s

plot(data(:, 1), data(:,2), 'r.')

axis( [-0.0, 1.1, 0, 1.1] )
xlabel('f_1')
ylabel('f_2')
title('MOEAD')
legend('ZDT1')
set(gca,'FontSize',16); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(get(gca,'XLabel'),'FontSize',16);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',16);
set(get(gca,'Title'),'FontSize',16);

%% ZDT2
clc,clear
clear all;
close all; 

data = load('ZDT2.txt');

plot(data(:, 1), data(:,2), 'r.')

axis( [0, 1.02, 0, 1.02] )
xlabel('f_1')
ylabel('f_2')
title('MOEAD')
legend('ZDT2')
set(gca,'FontSize',16); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(get(gca,'XLabel'),'FontSize',16);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',16);
set(get(gca,'Title'),'FontSize',16);

%% ZDT3
clc,clear
clear all;

data = load('ZDT3.txt');

plot(data(:, 1), data(:,2), 'r.')

axis( [0, 0.9, -0.9, 1.02] )
xlabel('f_1')
ylabel('f_2')
title('MOEAD')
legend('ZDT3')
set(gca,'FontSize',16); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(get(gca,'XLabel'),'FontSize',16);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',16);
set(get(gca,'Title'),'FontSize',16);

%% DYLZ1
clc,clear
clear all;

data = load('DTLZ1.txt');

plot3(data(:, 1), data(:,2),data(:,3) ,  'r.' )

axis( [0, 1, 0, 1, 0, 1] )
xlabel('f_1')
ylabel('f_2')
zlabel('f_3')
title('MOEAD')
legend('DTLZ1')
set(gca,'FontSize',16); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(get(gca,'XLabel'),'FontSize',16);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',16);
set(get(gca,'ZLabel'),'FontSize',16);
set(get(gca,'Title'),'FontSize',16);

%% DYLZ2
clc,clear
clear all;

data = load('DTLZ2.txt');

plot3(data(:, 1), data(:,2),data(:,3) ,  'r.' )

axis( [0, 1, 0, 1, 0, 1] )
xlabel('f_1')
ylabel('f_2')
zlabel('f_3')
title('MOEAD')
legend('DTLZ2')
set(gca,'FontSize',16); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(get(gca,'XLabel'),'FontSize',16);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',16);
set(get(gca,'ZLabel'),'FontSize',16);
set(get(gca,'Title'),'FontSize',16);

%% F1
clc,clear
clear all;

data = load('F4.txt'); 
% h = figure;
% set(h,'units','normalized','position',[0.2 0.2 0.5 0.5]);
figure(1)
plot(data(:, 1), data(:,2), 'r.')

axis( [-0.0, 1, 0, 1.1] )
xlabel('f_1')
ylabel('f_2')
title('MOEA/D-DE')
legend('F4-PF')
set(gca,'FontSize',16); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(get(gca,'XLabel'),'FontSize',16);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',16);
set(get(gca,'Title'),'FontSize',16);

data = load('F4_SOL.txt');

figure(2)
plot3(data(:, 1), data(:,2),data(:,3) ,  'r.' )

axis( [0,  1, -1, 1, -1, 1] )
xlabel('x1')
ylabel('x2')
zlabel('x3')
title('MOEA/D-DE')
legend('F4-SOL')
set(gca,'FontSize',16); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(get(gca,'XLabel'),'FontSize',16);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',16);
set(get(gca,'ZLabel'),'FontSize',16);
set(get(gca,'Title'),'FontSize',16);






