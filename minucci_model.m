% 
% clc
% clear all
%  close all
% 
% % ic_cal_exp = nan(1,39);
% % 
% %i.c. from Maiti model (latest times after running for 1000 hrs with no LPS)
% % ic_cal_exp(1) = 0;
% % ic_cal_exp(2) = Y0_maiti(end,5);
% % ic_cal_exp(3) = Y0_maiti(end,7);
% % ic_cal_exp(4) = Y0_maiti(end,15);
% % ic_cal_exp(5) = Y0_maiti(end,16);
% % ic_cal_exp(6) = Y0_maiti(end,20);
% % ic_cal_exp(7) = Y0_maiti(end,21);
% % ic_cal_exp(8) = Y0_maiti(end,17);
% % ic_cal_exp(9) = Y0_maiti(end,19);
% % ic_cal_exp(10) = Y0_maiti(end,18);
% % ic_cal_exp(11) = Y0_maiti(end,31);
% % ic_cal_exp(12) = Y0_maiti(end,11);
% % ic_cal_exp(13) = Y0_maiti(end,30);
% % ic_cal_exp(14) = Y0_maiti(end,12);
% % ic_cal_exp(15) = Y0_maiti(end,10);
% % ic_cal_exp(16) = Y0_maiti(end,6);
% % ic_cal_exp(17) = 1e-7;  %il10_rjt  
% % ic_cal_exp(18) = Y0_maiti(end,27);
% % ic_cal_exp(19) = Y0_maiti(end,22);
% % ic_cal_exp(20) = Y0_maiti(end,6);
% % ic_cal_exp(21) = 0; %jak1
% % ic_cal_exp(22) = Y0_maiti(end,3);
% % ic_cal_exp(23) = Y0_maiti(end,13);
% % ic_cal_exp(24) = Y0_maiti(end,14);
% % ic_cal_exp(25) = 0; %SOCS1cyto
% % ic_cal_exp(26) = 0;  %SOCS1mRNA
% % ic_cal_exp(27) = 0; %SOCS3cyto
% % ic_cal_exp(28) = 0;  %SOCS3mRNA
% % ic_cal_exp(29) = Y0_maiti(end,25); %STAT3a
% % ic_cal_exp(30) = Y0_maiti(end,24);  %STAT3i
% % ic_cal_exp(31) = Y0_maiti(end,26);
% % ic_cal_exp(32) = 1e-4;   %STAT3ni  
% % ic_cal_exp(33) = Y0_maiti(end,2);
% % ic_cal_exp(34) = Y0_maiti(end,9);
% % ic_cal_exp(35) = Y0_maiti(end,28);
% % ic_cal_exp(36) = Y0_maiti(end,23);
% % ic_cal_exp(37) = Y0_maiti(end,8);
% % ic_cal_exp(38) = 0;  %Tyk2
% % ic_cal_exp(39) = 0;  %il10 act
% % ic_cal_exp(39) = Y0_maiti(end,22) - Y0_maiti(end,14);
% 
% % %the ones included by them in the code
ic_cal_exp=[0.00E+00,1.25E-08,3.56E-09,5.90E-04,3.24E-08,2.53E-01,...
    5.18E-10,1.07E-01,1.58E-08,7.34E+00,4.24E-07,3.19E-10,2.69E-05,...
    1.47E-02,1.85E-01,2.33E-08,2.94E-10,2.85E-06,5.82E-08,1.00E-01,...
    2.00E-01,0.00E+00,3.97E-08,4.81E-09,2.59E-10,4.85E-08,5.86E-09,...
    1.38E-07,6.59E-09,3.50E-01,3.90E-09,4.20E-11,1.00E-01,4.40E-09,...
    1.40E-06,2.94E-08,1.00E-01,2.00E-01,5.20E-08]; 
% % 
params = [5.953e-4,8.736e-5,11.297,0.115,0.0714,3.57e-4,0.158,0.180,0.054,3.003e-5,0.0020,...
    5.719e-4,2.108e-4,0.275,0.04,0.0013,0.092,1.636e-4,0.0076,0.0246,0.335,0.229,0.007,...
    4.418e-4,0.938,0.0336,0.936,1.804e-5,0.0031,1.0192,1.97,0.0047,2.701,6.939e-5,1.427e-4,9.532e-5,...
    0.0081,0.0089,1.045,0.0134,0.0073,0.0234,3.5910,0.139,0.11,0.0717,0.009,0.0126,0.0364,0.237,...
    10.609,21.933,4.669e-6,0.388];
% 
% 
% %params included by them in the code
% 
% % params_cal_exp=[5.95E-05,8.74E-05,1.13E+01,1.15E-01,7.14E-02,3.57E-04,...
% %     1.58E-01,1.80E-01,1.54E-02,3.00E-05,2.02E-03,5.72E-04,2.11E-04,...
% %     2.75E-01,4.00E-02,1.34E-03,9.21E-02,1.64E-04,7.61E-03,2.46E-02,...
% %     3.42E-01,2.29E-01,7.05E-03,5.42E-04,9.38E-01,3.36E-02,9.36E-01,...
% %     1.80E-05,3.13E-03,1.92E-02,1.97E+00,4.67E-03,2.70E+00,6.94E-05,...
% %     1.43E-04,9.53E-05,8.08E-03,8.94E-03,1.05E+00,1.34E-02,7.27E-03,...
% %     2.34E-02,3.59E+00,1.39E-01,1.10E-01,7.17E-02,9.12E-03,1.26E-02,...
% %     3.64E-02,2.37E-01,1.06E+01,2.19E+01,4.67E-06,3.88E-01];
% 
% 
t_1 = 12*3600; %s
% 
ic_cal_exp(1)=10/75;  %LPS
% ic_cal_exp(3)=5.88e-4;  %TNF-a sup
% % ic_cal_exp(2)=5.56e-4;  %IL-10 sup
% 
[T_comp,Y_minucci_ode] = ode23s(@rhs_crosstalk_minucci,[0 t_1],ic_cal_exp,[],params,0);
% 
% 






%%%%%%% run ABM
abm_sim=abm_run;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot M1 & M2 activation - ABM & ODE | single macrophage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract ABM results from test
% a3: ODE and ABM modeling M1 and M2 activation
% extract transients from ABM called test

fig_name='cal exp';

tspan = 1:abm_sim.generations;     %da 1 al numero di iterazioni complessive  
tspan = tspan .* (abm_sim.GenerationSize/60);
tspan=[0 tspan];

avgm1act = abm_sim.Rule.Results{1}.avgm1act;
avgm2act = abm_sim.Rule.Results{1}.avgm2act;
avgfact = abm_sim.Rule.Results{1}.avgfact;

avgpimcount = abm_sim.Rule.Results{1}.avgpimcount;
avgaimcount = abm_sim.Rule.Results{1}.avgaimcount;
sdm1act = abm_sim.Rule.Results{1}.sdm1act;
sdm2act = abm_sim.Rule.Results{1}.sdm2act;
sdfact = abm_sim.Rule.Results{1}.sdfact;

sdpimcount = abm_sim.Rule.Results{1}.sdpimcount;
sdaimcount = abm_sim.Rule.Results{1}.sdaimcount;
avgm1count = abm_sim.Rule.Results{1}.avgm1act_tot;
sdm1countt = abm_sim.Rule.Results{1}.sdm1count;
avgm2count = abm_sim.Rule.Results{1}.avgm2act_tot;
sdm2countt = abm_sim.Rule.Results{1}.sdm2count;
avgintcount = abm_sim.Rule.Results{1}.avgintact_tot;
sdintcountt = abm_sim.Rule.Results{1}.sdmintcount;
avgm0count = abm_sim.Rule.Results{1}.avgm0count_tot;
sdm0countt = abm_sim.Rule.Results{1}.sdm0count;
avgfact_tot = abm_sim.Rule.Results{1}.avgfact_tot;
sdfcount = abm_sim.Rule.Results{1}.sdfcount;

dark_red=[168, 47, 75]/255;
dark_yellow=[176, 134, 35]/255;
dark_green=[3, 156, 115]/255;
dark_blue=[11, 93, 120]/255;
light_red=[255, 120, 151]/255;
light_yellow=[255, 209, 102]/255;
light_green=[6, 214, 159]/255;
light_blue=[42, 175, 219]/255;
% 
% % m1 activation
figure;
subplot(4,2,1)
% hold on
h1=boundedline(tspan,avgm1act,sdm1act,'alpha','cmap',dark_red,'linewidth',2);
% h2=plot(T_comp/3600,Y_comp(:,36)/max(Y_comp(:,36)),'linewidth',2,'color',dark_yellow);
ylabel('Average M1 activation')
xlabel('Time (h)')
xlim([0 24])
% legend([h1(1) h2(1)],'ABM','ODE')
set(gca,'fontsize',8)
% 
% % m2 activation
subplot(4,2,2)
% hold on
h1=boundedline(tspan,avgm2act,sdm2act,'alpha','cmap',light_red,'linewidth',2);
% h2=plot(T_comp/3600,Y_comp(:,19)/max(Y_comp(:,19)),'linewidth',2,'color',dark_green);
ylabel('Average M2 activation')
xlabel('Time (h)')
xlim([0 24])
% legend([h1(1) h2(1)],'ABM','ODE')
set(gca,'fontsize',8)
% 
% % pro-inflammatories/TNFa
subplot(4,2,7)
% hold on
h1=boundedline(tspan,avgpimcount,sdpimcount,'alpha','cmap',dark_blue,'linewidth',2); % light red
% h2=plot(T_comp/3600,(Y_comp(:,3))/max(Y_comp(:,3)),'color',light_yellow,'linewidth',2); % light yellow
% h3=plot(T_comp/3600,(Y_comp(:,3)+Y_comp(:,34))/max(Y_comp(:,3)+Y_comp(:,34)),'--','color',light_yellow,'linewidth',2);
ylabel('Average PIM count')
xlabel('Time (h)')
xlim([0 24])
% ylim([0 1])
% legend([h1(1) h2(1) h3(1)],'ABM','ODE','ODE (+ receptor-bound)')
set(gca,'fontsize',8)
% 
% % anti-inflammatories/IL-10
subplot(4,2,8)
% hold on
h1=boundedline(tspan,avgaimcount,sdaimcount,'alpha','cmap',light_blue,'linewidth',2); % light blue
% h2=plot(T_comp/3600,Y_comp(:,2)/max(Y_comp(:,2)),'color',light_green,'linewidth',2); % light green
% h3=plot(T_comp/3600,(Y_comp(:,2)+Y_comp(:,16))/max((Y_comp(:,2)+Y_comp(:,16))),'--','color',light_green,'linewidth',2); % light green
ylabel('Average AIM count')
xlabel('Time (h)')
xlim([0 24])
% ylim([0 1])
% legend([h1(1) h2(1) h3(1)],'ABM','ODE','ODE (+ receptor-bound)')
set(gca,'fontsize',8)

subplot(4,2,3)
h1=boundedline(tspan,avgm1count,sdm1countt,'alpha','cmap',dark_yellow,'linewidth',2); % light blue
ylabel('Average M1 count')
xlabel('Time (h)')
xlim([0 24])
% ylim([0 1])
% legend([h1(1) h2(1) h3(1)],'ABM','ODE','ODE (+ receptor-bound)')
set(gca,'fontsize',8)

subplot(4,2,4)
h1=boundedline(tspan,avgm2count,sdm2countt,'alpha','cmap',light_yellow,'linewidth',2); % light blue
ylabel('Average M2 count')
xlabel('Time (h)')
xlim([0 24])
% ylim([0 1])
% legend([h1(1) h2(1) h3(1)],'ABM','ODE','ODE (+ receptor-bound)')
set(gca,'fontsize',8)


subplot(4,2,5)
% hold on
h1=boundedline(tspan,avgfact_tot,sdfcount,'alpha','cmap',light_green,'linewidth',2); % light red
% h2=plot(T_comp/3600,(Y_comp(:,3))/max(Y_comp(:,3)),'color',light_yellow,'linewidth',2); % light yellow
% h3=plot(T_comp/3600,(Y_comp(:,3)+Y_comp(:,34))/max(Y_comp(:,3)+Y_comp(:,34)),'--','color',light_yellow,'linewidth',2);
ylabel('Average F count')
xlabel('Time (h)')
xlim([0 24])
% ylim([0 1])
% legend([h1(1) h2(1) h3(1)],'ABM','ODE','ODE (+ receptor-bound)')
set(gca,'fontsize',8)
% 
% % anti-inflammatories/IL-10
subplot(4,2,6)
% hold on
h1=boundedline(tspan,avgfact,sdfact,'alpha','cmap',dark_green,'linewidth',2); % light blue
% h2=plot(T_comp/3600,Y_comp(:,2)/max(Y_comp(:,2)),'color',light_green,'linewidth',2); % light green
% h3=plot(T_comp/3600,(Y_comp(:,2)+Y_comp(:,16))/max((Y_comp(:,2)+Y_comp(:,16))),'--','color',light_green,'linewidth',2); % light green
ylabel('Average F act')
xlabel('Time (h)')
xlim([0 24])
% ylim([0 1])
% legend([h1(1) h2(1) h3(1)],'ABM','ODE','ODE (+ receptor-bound)')
set(gca,'fontsize',8)
% 
% figure('Name',fig_name)
% hold on
% h1=boundedline(tspan,avgm1act/max(avgm1act),sdm1act/max(avgm1act),'alpha','cmap',dark_red,'linewidth',2);
% h2=plot(T_comp,Y_comp(:,36)/max(Y_comp(:,36)),'linewidth',2,'color',dark_yellow);
% h3=boundedline(tspan,avgm2act/max(avgm1act),sdm2act/max(avgm1act),'alpha','cmap',dark_blue,'linewidth',2);
% h4=plot(T_comp,Y_comp(:,19)/max(Y_comp(:,36)),'linewidth',2,'color',dark_green);
% ylabel('Average M1/M2 activation')
% xlabel('hours')
% xlim([0 t_1])
% ylim([0 1])
% leg1([h1(1) h2(1) h3(1) h4(1)],'M1 (ABM)','M1 (ODE)','M2 (ABM)','M2 (ODE)')
% set(gca,'fontsize',8)
% 
% 
% 
% 
% 
% 
