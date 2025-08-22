clc
clear all 
close all
% import_Data =xlsread('DataSet-AliAkbar.xlsx','B3:C6274')
  import_Data1=xlsread('F:\رسالة ماجستير\عملي\MATLA\codes\DataSet-AliAkbar',2,'B3:C6274')
    Data.MSE_encrypted =import_Data1(:,2);
    Data.MSE_decrypted =import_Data1(:,3);
   % Data.V=import_data(:,5);
    %Data.SOC=import_data(:,6);

    %% Plot Data
    figure
    plot( Data.MSE_encrypted, Data.MSE_decrypted)
  %  figure
   % plot(Data.T,Data.V)
    %figure