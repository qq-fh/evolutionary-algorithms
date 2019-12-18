%% ����1
clc,clear
close all

fitness = load('optiFitnesss_func1.txt');
plot(fitness, 'r')
axis( [-3000, 51000, -13000, -1500] )
xlabel('Number of Generation')
ylabel('Value')
legend('func1')
set(gca,'FontSize',16); % �������ִ�С��ͬʱӰ���������ע��ͼ��������ȡ�
set(get(gca,'XLabel'),'FontSize',16);%ͼ������Ϊ8 point��С5��
set(get(gca,'YLabel'),'FontSize',16);



%% ����3
clc,clear
close all

fitness = load('optiFitnesss_func3.txt');
plot(fitness, 'r')
axis( [-1000, 21000, -1, 22] )
xlabel('Number of Generation')
ylabel('Value')
legend('func3')
set(gca,'FontSize',16); % �������ִ�С��ͬʱӰ���������ע��ͼ��������ȡ�
set(get(gca,'XLabel'),'FontSize',16);%ͼ������Ϊ8 point��С5��
set(get(gca,'YLabel'),'FontSize',16);


%% ����4
clc,clear
close all

fitness = load('optiFitnesss_func4.txt');
plot(fitness, 'r')
axis( [-100, 1100, -50, 700] )
xlabel('Number of Generation')
ylabel('Value')
legend('func4')
set(gca,'FontSize',16); % �������ִ�С��ͬʱӰ���������ע��ͼ��������ȡ�
set(get(gca,'XLabel'),'FontSize',16);%ͼ������Ϊ8 point��С5��
set(get(gca,'YLabel'),'FontSize',16);

%% ����9
clc,clear
close all

fitness = load('optiFitnesss_func9.txt');
plot(fitness, 'r')
axis( [-5, 105, -10, 120] )
xlabel('Number of Generation')
ylabel('Value')
legend('func4')
set(gca,'FontSize',16); % �������ִ�С��ͬʱӰ���������ע��ͼ��������ȡ�
set(get(gca,'XLabel'),'FontSize',16);%ͼ������Ϊ8 point��С5��
set(get(gca,'YLabel'),'FontSize',16);

