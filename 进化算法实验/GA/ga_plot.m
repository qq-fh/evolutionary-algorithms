%% 函数1
clc,clear
close all

fitness = load('optiFitnesss_func1.txt');
plot(fitness, 'r')
axis( [-3000, 51000, -13000, -1500] )
xlabel('Number of Generation')
ylabel('Value')
legend('func1')
set(gca,'FontSize',16); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(get(gca,'XLabel'),'FontSize',16);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',16);



%% 函数3
clc,clear
close all

fitness = load('optiFitnesss_func3.txt');
plot(fitness, 'r')
axis( [-1000, 21000, -1, 22] )
xlabel('Number of Generation')
ylabel('Value')
legend('func3')
set(gca,'FontSize',16); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(get(gca,'XLabel'),'FontSize',16);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',16);


%% 函数4
clc,clear
close all

fitness = load('optiFitnesss_func4.txt');
plot(fitness, 'r')
axis( [-100, 1100, -50, 700] )
xlabel('Number of Generation')
ylabel('Value')
legend('func4')
set(gca,'FontSize',16); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(get(gca,'XLabel'),'FontSize',16);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',16);

%% 函数9
clc,clear
close all

fitness = load('optiFitnesss_func9.txt');
plot(fitness, 'r')
axis( [-5, 105, -10, 120] )
xlabel('Number of Generation')
ylabel('Value')
legend('func4')
set(gca,'FontSize',16); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(get(gca,'XLabel'),'FontSize',16);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',16);

