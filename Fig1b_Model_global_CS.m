close all
clear all

addpath('./Utility/')

cell_subtype = {'PC', 'PV', 'SST', 'VIP'};
color_subtype = {[1,0,0], [0,0,1], [1,.5,.0], [.7,.7,.7]};
Nsub = length(cell_subtype);

% - population size
N_pc = 500;
N_pv = 100;
N_st = 100;
N_vip= 100;

ids_pc = 1:N_pc;
ids_pv = N_pc+1:N_pc+N_pv;
ids_st = N_pc+N_pv+1:N_pc+N_pv+N_st;
ids_vip = N_pc+N_pv+N_st+1:N_pc+N_pv+N_st+N_vip;

ids_all{1} = ids_pc;
ids_all{2} = ids_pv;
ids_all{3} = ids_st;
ids_all{4} = ids_vip;

N = N_pc + N_pv + N_st +N_vip;

% - spatial distribution
L = 1; %mm

x_pc = 0:L/N_pc:L-L/N_pc;
x_pv = 0:L/N_pv:L-L/N_pv;
x_st = 0:L/N_st:L-L/N_st;
x_vip = 0:L/N_vip:L-L/N_vip;

x_all{1} = x_pc;
x_all{2} = x_pv;
x_all{3} = x_st;
x_all{4} = x_vip;

% - simulation parameters (dynamics + plasticity)
dt = 1;
tau = 10.;

th = zeros(1,N);
%th(ids_pc) = 0.5;
%th(ids_pv) = 0.5;

t_sim = 1000;
T = 0:dt:t_sim;
t_init = 100;
t0 = 200;
t1 = 400;
t2 = 600;
t3 = 800;
t4 = 1000;
ttr=50;

% - Baseline input, signal amplitude for [PN,PV,SST,VIP]
amp= 2;                             %amplitude
bb = [0,0,0,0];                     %baseline 
mm = [amp, 2*amp, 2*amp, 2*amp];    %max amplitude

% - Spatial tuning setting
alpha = [1,1,1,1];                        

%   - Input heterogeneity (toggle one)
%1) - no heterogenity
rl= 1;
ru= 0;

% %2) - with heterogenity
% rl= 0
% ru= 2

%% - weight matrix: CS effects

% - global connectivity strength: CS
A = [.1,.33,.33, .33;               %PN->  PN PV SST VIP
    -.3,-.57,-.28, -0.08;            %PV->  PN PV SST VIP
    -.56, -.35, -.29, -0.37;        %SST-> PN PV SST VIP
    -0.04, -.13, -.15, -.06];      %VIP-> PN PV SST VIP

% - spatial connectivity parameters: CS
S = [0.15, 0.15, 0.15, 0.15; 
     0.13, 0.15, 0.14, 0.22; 
     0.21, 0.17, 0.33, 0.2; 
     0.19, 0.14, 0.21, 0.27]; 

% - input center offset (according to medio-lateral connectivity maps): CS
 d0_matrix = [0, 0, 0, 0                     
             0.001, 0.032, -0.005, 0.026    
             0.037, 0.024, 0.079, 0.038     
             0.013, -0.033, 0.016, 0.068];         

Wi = [];
Di = [];
for i = 1:Nsub
    for j = 1:Nsub
        d0 = d0_matrix(i,j);

        %[Wi{i,j}, Di{i,j}] = dd_w_fun_2(x_all{i}, x_all{j}, L, A(i,j), S(i,j));
        [Wi{i,j}, Di{i,j}] = dd_w_fun(x_all{i}, x_all{j}, L, A(i,j), S(i,j), d0);
    end
end

%% This section tests the inhibitory effects of PV and SST and disinhibitory effects of VIP connections


% % - IN motif _________________  CHECK INTERNEURON_________________________
%fig_suffix = ['global CS PV_SST inh and VIP dis_ alpha_', num2str(alpha(1))]  
fig_suffix = 'Fig1b_Model_global_CS';

% %inhibition
W1 = [Wi{1,1}*1, Wi{1,2}*0, Wi{1,3}*0, Wi{1,4}*0;     
     Wi{2,1}*1, Wi{2,2}*1, Wi{2,3}*1, Wi{2,4}*1;      
     Wi{3,1}*1, Wi{3,2}*1, Wi{3,3}*1, Wi{3,4}*1;      
     Wi{4,1}*0, Wi{4,2}*1, Wi{4,3}*1, Wi{4,4}*0;      
     ]; 

W2 = [Wi{1,1}*1, Wi{1,2}*0, Wi{1,3}*0, Wi{1,4}*0; 
     Wi{2,1}*0, Wi{2,2}*0, Wi{2,3}*0, Wi{2,4}*0;  
     Wi{3,1}*0, Wi{3,2}*0, Wi{3,3}*0, Wi{3,4}*0;   
     Wi{4,1}*0, Wi{4,2}*0, Wi{4,3}*0, Wi{4,4}*0; 
      ]; 

%% - inputs
I0 = .1+.25/10*(rand([N,length(T)])-.5);
I0 = I0*0;

x_stim = L/2;

Im = zeros(N,length(T));
for k = 1:Nsub
    for kk = 1:length(ids_all{k})
        d = x_all{k}(kk) - x_stim;
        
        if k == 1
        Im(ids_all{k}(kk),(t0:t3)/dt) = bb(k)+mm(k)*exp(alpha(k)*(cos(d/L*2*pi)))/exp(alpha(k))*(rl + rand() * ru);
        end

        if k == 2
        Im(ids_all{k}(kk),(t1:t3)/dt) = bb(k)+mm(k)*exp(alpha(k)*(cos(d/L*2*pi)))/exp(alpha(k))*(rl + rand() * ru);
        end

        if k == 3
        Im(ids_all{k}(kk),(t1:t3)/dt) = bb(k)+mm(k)*exp(alpha(k)*(cos(d/L*2*pi)))/exp(alpha(k))*(rl + rand() * ru);
        end

        if k == 4
        Im(ids_all{k}(kk),(t2:t3)/dt) = bb(k)+mm(k)*exp(alpha(k)*(cos(d/L*2*pi)))/exp(alpha(k))*(rl + rand() * ru);
        end

    end
end

I = rectify(I0 + Im); 

%% - response dynamics for both network simulations 

[r1,inp1] = simulate_network2(I, W1, T, dt, tau, th);

[r2,inp2] = simulate_network2(I, W2, T, dt, tau, th);


%% - figures

ts1 = (t0+t1)/2;
ts2 = (t1+t2)/2;
ts3 = (t2+t3)/2;
ts4 = (t3+t4)/2;

figure()

subplot(121)
hold on
title(['t: ' num2str(ts2/dt)]);
% r0 = (r_all(ids_pc,ts2/dt) + r_all(ids_pc,ts4/dt))/2;

plot(r1(ids_pc,ts2/dt),'linewidth', 1)
plot(r2(ids_pc,ts2/dt),'linewidth', 1)

plot(r1(ids_pc,ts2/dt) - r2(ids_pc,ts2/dt),'linewidth', 1)
activityIN_n=r1(ids_pc,ts2/dt);                                               %neuronal activity
activity_baseline_n=r2(ids_pc,ts2/dt);
act_difference_n=r1(ids_pc,ts2/dt) - r2(ids_pc,ts2/dt);

xlabel('Neuron #')
ylabel('Activity')

legend({'r1','r2', 'r1-r2'}); %['t: ' num2str(ts1)], ['t: ' num2str(ts3)], ['t: ' num2str(ts5)]}, 'Location','best')
legend boxoff

set(gca, 'LineWidth', 1, 'FontSize', 15, 'TickDir', 'out', 'TickLength',[.025,.025]/2)

subplot(122)
hold on
title(['t: ' num2str(ts2/dt)]);

plot(inp1(ids_pc,ts2/dt),'linewidth', 1)
plot(inp2(ids_pc,ts2/dt),'linewidth', 1)

plot(inp1(ids_pc,ts2/dt) - inp2(ids_pc,ts2/dt),'linewidth', 1)
activityIN=inp1(ids_pc,ts2/dt);
activity_baseline=inp2(ids_pc,ts2/dt);
act_difference=inp1(ids_pc,ts2/dt) - inp2(ids_pc,ts2/dt);

xlabel('Neuron #')
ylabel('Activity')

legend({'inp 1','inp 2', 'inp 1 - inp 2'}); %['t: ' num2str(ts1)], ['t: ' num2str(ts3)], ['t: ' num2str(ts5)]}, 'Location','best')
legend boxoff

set(gca, 'LineWidth', 1, 'FontSize', 15, 'TickDir', 'out', 'TickLength',[.025,.025]/2)

subfolderName = [fig_suffix]; % Modify this line

% Check if the subfolder exists
if ~exist(subfolderName, 'dir')
    % If the subfolder doesn't exist, create it
    mkdir(subfolderName);
end

% Define the file name within the "output" subfolder
%fig_suffix = 'my_figure';  % Replace with your desired file name or suffix
outputFilePath = fullfile(subfolderName, ['/a.png']);

% Save the figure in the "output" subfolder
print(outputFilePath, '-dpng', '-r300');


%% - everything in one
plot_everything = 1;
r_all = inp1;

if plot_everything

figure('Position',[100,100,1200,600])

subplot(231); hold on; title('Weight matrix')
imagesc(W1)
colorbar()
xlabel('Neuron #')
ylabel('Neuron #')
axis image

set(gca, 'LineWidth', 1, 'FontSize', 15, 'TickDir', 'out', 'TickLength',[.025,.025]/2)

subplot(232); hold on; title('Input')
imagesc(I, [0,4.5])
colorbar()
xlim([0,t_sim/dt])
ylim([1,N])

xlabel('Time')
ylabel('Neuron #')

set(gca, 'LineWidth', 1, 'FontSize', 15, 'TickDir', 'out', 'TickLength',[.025,.025]/2)

subplot(233); hold on; title('Activity')
imagesc(r_all, [0,4.5])
colorbar()
xlim([0,t_sim/dt])
ylim([1,N])

xlabel('Time')
ylabel('Neuron #')

set(gca, 'LineWidth', 1, 'FontSize', 15, 'TickDir', 'out', 'TickLength',[.025,.025]/2)

subplot(234); hold on; %title('Activity')
%plot(r_all(1,:), 'Color',color_subtype{1})
plot(T, r_all(N_pc/2,:), 'linewidth', 1, 'Color',color_subtype{1})

%plot(r_all(N_pc+1,:), 'Color',color_subtype{2})
plot(T, r_all(N_pc+N_pv/2,:), 'linewidth', 1,'Color',color_subtype{2})

%plot(r_all(N_pc+N_pv+1,:), 'Color',color_subtype{3})
plot(T, r_all(N_pc+N_pv+N_st/2,:), 'linewidth', 1,'Color',color_subtype{3})

plot(T, r_all(N_pc+N_pv+N_st+N_vip/2,:), 'linewidth', 1,'Color',color_subtype{4})

legend(cell_subtype, 'Location','best');%, 'NumColumns',3)
legend boxoff

xlim([0,t_sim])
xlabel('Time')
ylabel('Activity')

%ylim([0,2])

set(gca, 'LineWidth', 1, 'FontSize', 15, 'TickDir', 'out', 'TickLength',[.025,.025]/2)

subplot(235); hold on; %title('Activity')

% plot(r_all(:,ts1/dt),'linewidth', 1)
plot(r_all(:,ts2/dt),'linewidth', 1)
% plot(r_all(:,ts3/dt),'linewidth', 1)
% plot(r_all(:,ts4/dt),'linewidth', 1)


legend({['t: ' num2str(ts1)], ['t: ' num2str(ts2)], ['t: ' num2str(ts3)], ['t: ' num2str(ts4)]}, 'Location','best')
legend boxoff

xlim([0,N])

xlabel('Neuron #')
%ylabel('Activity')

%ylim([0,2])

set(gca, 'LineWidth', 1, 'FontSize', 15, 'TickDir', 'out', 'TickLength',[.025,.025]/2)

subplot(236); hold on; %title('Activity')

plot(I(:,ts1/dt))
plot(I(:,ts2/dt))
plot(I(:,ts3/dt))
plot(I(:,ts4/dt))

legend({['t: ' num2str(ts1)], ['t: ' num2str(ts2)], ['t: ' num2str(ts3)], ['t: ' num2str(ts4)]}, 'Location','best')
legend boxoff

xlim([0,N])

xlabel('Neuron #')
ylabel('Input')

%ylim([0,2])

set(gca, 'LineWidth', 1, 'FontSize', 15, 'TickDir', 'out', 'TickLength',[.025,.025]/2)

%print(['./Allin1_' fig_suffix '.png'], '-dpng', '-r300');

% Define the subfolder name with fig_suffix
subfolderName = [fig_suffix]; % Modify this line

% Check if the subfolder exists
if ~exist(subfolderName, 'dir')
    % If the subfolder doesn't exist, create it
    mkdir(subfolderName);
end

% Define the file name within the "output" subfolder
%fig_suffix = 'my_figure';  % Replace with your desired file name or suffix
outputFilePath = fullfile(subfolderName, ['./b.png']);

% Save the figure in the "output" subfolder
print(outputFilePath, '-dpng', '-r300');
end

