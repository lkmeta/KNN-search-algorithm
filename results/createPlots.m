function [] = createPlots(filename)
%createPlots.m Function creating the plots shown in the report and the
%diagrams folder on github.
%Input:
%   filename: the name of the matrix
%Output:
%   None

%create the file path
filepath0 = strcat(filename,'v0.txt');
filepath1 = strcat(filename,'v1.txt');
filepath2 = strcat(filename,'v2.txt');

%store string with n,d
if strcmp(filename,'BBC')
    sizestr = ' n=17720, d=17';
    
elseif strcmp(filename,'CNN')
    sizestr = ' n=22545, d=17';
    
elseif strcmp(filename,'CNNIBN')
    sizestr = ' n=33117, d=17';
    
elseif strcmp(filename,'NDTV')
    sizestr = ' n=17051, d=17';
    
elseif strcmp(filename,'TIMESNOW')
    sizestr = ' n=39252, d=17';
    
elseif strcmp(filename,'MiniBooNE_PID')
    sizestr = ' n=130064, d=50';
    
elseif strcmp(filename,'ColorHistogram')
    sizestr = ' n=68040, d=32';
    
elseif strcmp(filename,'ColorMoments')
    sizestr = ' n=68040, d=9';
    
elseif strcmp(filename,'CoocTexture')
    sizestr = ' n=68040, d=16';
    
elseif strcmp(filename,'LayoutHistogram')
    sizestr = ' n=62481, d=32';
    
end

%import data
data0 = importdata(filepath0);
data1 = importdata(filepath1);
data2 = importdata(filepath2);


%v0: time - k
x = data0(1,2:5);
y = data0(2,2:5);

titlestr = strcat('v0 for',{' '},filename,':',sizestr);

figure
plot(x,y,'-*b')
xlabel('k','fontweight','bold')
ylabel('execution time(s)','fontweight','bold')
title(titlestr)

%v1: time - processes
x = data1(2:6,1);
y = data1(2:6,2:5);

titlestr = strcat('v1 for',{' '},filename,':',sizestr);

figure
plot(x,y(:,1),'-b')
hold on
plot(x,y(:,2),'-r')
plot(x,y(:,3),'-g')
plot(x,y(:,4),'-y')
hold off
xlabel('# of processes','fontweight','bold')
ylabel('execution time(s)','fontweight','bold')
title(titlestr)
legend('k=10','k=20','k=50','k=100')

%v1: acceleration - processes
x = data1(2:6,1);
y = data1(2:6,2:5);
y0 = data0(2,2:5);

titlestr = strcat('v1 for',{' '},filename,':',sizestr);

for j=1:4
    for i=1:5
        y(i,j) = y0(j)/y(i,j);
    end
end

figure
plot(x,y(:,1),'-b')
hold on
plot(x,y(:,2),'-r')
plot(x,y(:,3),'-g')
plot(x,y(:,4),'-y')
hold off
xlabel('# of processes','fontweight','bold')
ylabel('acceleration','fontweight','bold')
title(titlestr)
legend('k=10','k=20','k=50','k=100')

%v1: time - k
x = data1(1,2:5);
y1 = data1(2,2:5);
y2 = data1(6,2:5);

titlestr = strcat('v1 for',{' '},filename,':',sizestr);

figure
plot(x,y1,'-b')
hold on
plot(x,y2,'-r')
hold off
xlabel('k','fontweight','bold')
ylabel('execution time(s)','fontweight','bold')
title(titlestr)
legend('processes=2','processes=20')

%v2: time - processes
x = data2(2:6,1);
y = data2(2:6,2:5);

titlestr = strcat('v2 for',{' '},filename,':',sizestr);

figure
plot(x,y(:,1),'-b')
hold on
plot(x,y(:,2),'-r')
plot(x,y(:,3),'-g')
plot(x,y(:,4),'-y')
hold off
xlabel('# of processes','fontweight','bold')
ylabel('execution time(s)','fontweight','bold')
title(titlestr)
legend('k=10','k=20','k=50','k=100')

%v2: acceleration - processes
x = data2(2:6,1);
y = data2(2:6,2:5);
y0 = data0(2,2:5);

titlestr = strcat('v2 for',{' '},filename,':',sizestr);

for j=1:4
    for i=1:5
        y(i,j) = y0(j)/y(i,j);
    end
end

figure
plot(x,y(:,1),'-b')
hold on
plot(x,y(:,2),'-r')
plot(x,y(:,3),'-g')
plot(x,y(:,4),'-y')
hold off
xlabel('# of processes','fontweight','bold')
ylabel('acceleration','fontweight','bold')
title(titlestr)
legend('k=10','k=20','k=50','k=100')

%v2: time - k
x = data2(1,2:5);
y1 = data2(2,2:5);
y2 = data2(6,2:5);

titlestr = strcat('v2 for',{' '},filename,':',sizestr);

figure
plot(x,y1,'-b')
hold on
plot(x,y2,'-r')
hold off
xlabel('k','fontweight','bold')
ylabel('execution time(s)','fontweight','bold')
title(titlestr)
legend('processes=2','processes=20')

%v1 vs v2 times - processes (k=50)
x = data1(2:6,1);
y1 = data1(2:6,4);
y2 = data2(2:6,4);

titlestr = strcat('v1 vs v2 for',{' '},filename,':',sizestr);

figure
plot(x,y1,'-b')
hold on
plot(x,y2,'-r')
hold off
xlabel('# of processes','fontweight','bold')
ylabel('execution time(s)','fontweight','bold')
title(titlestr)
legend('v1','v2')

end
