close all; clc; clearvars -except Data3

files = dir('*.txt');
for i=1:length(files)
    eval(['load ' files(i).name ' -ascii']);
end

%this line of code puts all the tables in a seperate matrix. To do this you
%have to copy the file names of the datasets from the workspace (right handed column) into this
%table

Voltage_exp= [V301, V351, V401, V451, V501, V551, V601, V651, V701, V751, V801, V851, V901, V951, V1001, V1051, V1101, V1151, V1201, V1251, V1301, V1351, V1401, V1451, V1501, V1551, V1601, V1651, V1701, V1751, V1801, V1851, V1901, V1951, V2001, V2051, V2101, V2151, V2201, V2251, V2301, V2351, V1401];
Voltage_exp2 = [V302, V352, V402, V452, V502, V552, V602, V652, V702, V752, V802, V852, V902, V952, V1002, V1052, V1102, V1152, V1202, V1252, V1302, V1352, V1402, V1452, V1502, V1552, V1602, V1652, V1702, V1752, V1802, V1852, V1902, V1952, V2002, V2052, V2102, V2152, V2202, V2252, V2302, V2352, V1402];
Voltage_exp3 = [V30, V35, V40, V45, V50, V55, V60, V65, V70, V75, V80, V85, V90, V95, V100, V105, V110, V115, V120, V125, V130, V135, V140, V145, V150, V155, V160, V165, V170, V175, V180, V185, V190, V195, V200, V205, V210, V215, V220, V225, V230, V235, V140];

%% for loop dataset 1
Sampling_Rate = 2500;                      
Measuring_Time = 1;
Length_Signal = Sampling_Rate * Measuring_Time;
n=1;
for i=2:2:size(Voltage_exp,2)
Data_raw2(:,n) = [Voltage_exp(:,i)];                                          %only consider the non-sample columns and put it in a table

Mean(:,n)  = mean(Data_raw2(:,n));                                            %calculate mean of signal in order to compensate for translation on y-axis
Data_02(:,n) = Data_raw2(:,n) - Mean(:,n);                                    %calculate the difference between mean and the data
B2(:,n) = fft(Data_02(:,n));                                                  %FFT to determine how many points we need to consider for the maximum based on frequency of signal
max_f(:,n)= max(abs(B2(:,n)));                                                %find the frecuency of the system at every voltage
amount_peaks(:,n)= uint64(max_f(:,n).*Sampling_Rate*10^-3*2);                 %calculate the integer of the ammount of peaks 
%integer_peaks(:,n)= uint64(ammount_peaks(:,n));
sort_data(:,n)= sort(abs(Data_02(:,n)),'desc');                               %sort the data so the top ammount of peaks are listed at the top of the table
high_peaks{n}(:,:)= sort_data(1:amount_peaks(:,n),n);                         %only take the highest peaks to analyse
mean_high_peaks(:,n)= mean(high_peaks{n}(:,:));                               %calcultate the mean of the highest peaks

Intensities_exp(:,n) = (mean_high_peaks(:,n).*104./10); %.*10e-3

n=n+1;
end

%% for loop dataset 2
clear Sampling_Rate Measuring_Time Length_Signal n i Data_raw2 Mean Data_02 B2 max_f amount_peaks sort_data high_peaks mean_high_peaks;
Sampling_Rate = 2500;                      
Measuring_Time = 1;
Length_Signal = Sampling_Rate * Measuring_Time;
n=1;
for i=2:2:size(Voltage_exp2,2)
Data_raw2(:,n) = [Voltage_exp2(:,i)];                                         %only consider the non-sample columns and put it in a table

Mean(:,n)  = mean(Data_raw2(:,n));                                            %calculate mean of signal in order to compensate for translation on y-axis
Data_02(:,n) = Data_raw2(:,n) - Mean(:,n);                                    %calculate the difference between mean and the data
B2(:,n) = fft(Data_02(:,n));                                                  %FFT to determine how many points we need to consider for the maximum based on frequency of signal
max_f(:,n)= max(abs(B2(:,n)));                                                %find the frecuency of the system at every voltage
amount_peaks(:,n)= uint64(max_f(:,n).*Sampling_Rate*10^-3*2);                 %calculate the integer of the ammount of peaks 
%integer_peaks(:,n)= uint64(ammount_peaks(:,n));
sort_data(:,n)= sort(abs(Data_02(:,n)),'desc');                               %sort the data so the top ammount of peaks are listed at the top of the table
high_peaks{n}(:,:)= sort_data(1:amount_peaks(:,n),n);                         %only take the highest peaks to analyse
mean_high_peaks(:,n)= mean(high_peaks{n}(:,:));                               %calcultate the mean of the highest peaks

Intensities_exp2(:,n) = (mean_high_peaks(:,n).*104./10); %.*10e-3

n=n+1;
end

%% for loop dataset 3
clear Sampling_Rate Measuring_Time Length_Signal n i Data_raw2 Mean Data_02 B2 max_f amount_peaks sort_data high_peaks mean_high_peaks;
Sampling_Rate = 2500;                      
Measuring_Time = 1;
Length_Signal = Sampling_Rate * Measuring_Time;
n=1;
for i=2:2:size(Voltage_exp3,2)
Data_raw2(:,n) = [Voltage_exp3(:,i)];                                         %only consider the non-sample columns and put it in a table

Mean(:,n)  = mean(Data_raw2(:,n));                                            %calculate mean of signal in order to compensate for translation on y-axis
Data_02(:,n) = Data_raw2(:,n) - Mean(:,n);                                    %calculate the difference between mean and the data
B2(:,n) = fft(Data_02(:,n));                                                  %FFT to determine how many points we need to consider for the maximum based on frequency of signal
max_f(:,n)= max(abs(B2(:,n)));                                                %find the frecuency of the system at every voltage
amount_peaks(:,n)= uint64(max_f(:,n).*Sampling_Rate*10^-3*2);                 %calculate the integer of the ammount of peaks 
%integer_peaks(:,n)= uint64(ammount_peaks(:,n));
sort_data(:,n)= sort(abs(Data_02(:,n)),'desc');                               %sort the data so the top ammount of peaks are listed at the top of the table
high_peaks{n}(:,:)= sort_data(1:amount_peaks(:,n),n);                         %only take the highest peaks to analyse
mean_high_peaks(:,n)= mean(high_peaks{n}(:,:));                               %calcultate the mean of the highest peaks

Intensities_exp3(:,n) = (mean_high_peaks(:,n).*104./10); %.*10e-3

n=n+1;
end


%% data uit het model
disp("if there's a error direct below, import Data3, the data from the model, in Matlab and change data2 -> data3 in the import venster")
Intensities_mod = table2array(Data3);
Intensities_mod(:,2) = Intensities_mod(:,2)*1000;
disp("Data from model is successfully imported")

%Intensities_mod2(:,1) = Intensities_mod(:,1);
%Intensities_mod2(:,2) = Intensities_mod(:,2)*1000;
%Intensities_mod = Intensities_mod2
%% plot alles in een
Voltage_plot = 3.0:0.5:24.0;
figure(1)
plot(Voltage_plot,Intensities_exp)
hold on
plot(Voltage_plot,Intensities_exp2)
plot(Voltage_plot,Intensities_exp3)
plot(Voltage_plot,Intensities_mod(:,2))
ylabel('Mean amplitude [mm]')
xlabel('Voltage[V]')
legend("experiment 1","experiment 2","experiment 3","matlab model") 

%% plot verschil tussen experiments en 
n = 30;
for i=1:1:43
    difference1(i,:) = abs(Intensities_exp(:,i) - Intensities_mod(i,2));
    difference2(i,:) = abs(Intensities_exp2(:,i) - Intensities_mod(i,2));
    difference3(i,:) = abs(Intensities_exp3(:,i) - Intensities_mod(i,2));
    
    %difference1(i,:) = Intensities_mod(i,2);%-Intensities_exp(:,i)
    %difference2(i,:) = Intensities_mod(i,2);%-Intensities_exp2(:,i) 
    %difference3(i,:) = Intensities_mod(i,2);%-Intensities_exp3(:,i)
    
    n = n+5;
end
figure(2)
plot(Voltage_plot,difference1)
hold on
plot(Voltage_plot,difference2)
plot(Voltage_plot,difference3)
ylabel('Difference between in [mm]')
xlabel('Voltage[V]')
legend("experiment 1","experiment 2","experiment 3") 

%% %error = (accepted - model)/ accepted * 100
clear errors
for i = 1:1:43;
    % https://www.sophia.org/tutorials/accuracy-and-precision--3
    
    error1 = difference1(i,:) / Intensities_exp(:,i) * 100;
    error2 = difference2(i,:) / Intensities_exp2(:,i) * 100;
    error3 = difference3(i,:) / Intensities_exp3(:,i) * 100;
    
    error_total1(i,1) = error1;
    error_total2(i,1) = error2;
    error_total3(i,1) = error3;
end

figure(3)
plot(Voltage_plot,error_total1)
hold on
plot(Voltage_plot,error_total2)
plot(Voltage_plot,error_total3)
ylabel('Error [%]')
xlabel('Voltage[V]')
legend("experiment 1","experiment 2","experiment 3") 

mean_error1 = mean(error_total1);
mean_error2 = mean(error_total2);
mean_error3 = mean(error_total3);

mean_total = (mean_error1+mean_error2+mean_error3)/3;

accuracy = 100 - mean_total

disp("Validation is done, there have to appear 3 graphs")
