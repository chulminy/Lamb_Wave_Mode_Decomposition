%% Validation of the proposed technique using FEM simulation
%
%% Contact
%  Name  : Chul Min Yeum
%  Email : chulminy@gmail.com
%  Please contact me if you have a question or find a bug

%% Simulation configuration
% 
% Please refer a detailed description of FEM simulation in Section 3.1 A
% major different is the dual PZT configuration on the plate. In the paper,
% dual PZTs A and B are collocated so as to generate individual S0 and A0
% modes. On the other hand, the FEM simulation used in this code is dual
% PZTs B and C are collocated to decompose S0 and A0 signla. Regardless of
% changes in dual PZT configuration, validation of the method is still
% valide due to reciprocity. 
%
%  Here is locations of dual PZTs A, B, and C
%  Dual PZT A            : Excitation
%  Dual PZTs B and C     : Sensing 
% 
%   A                   B
% -------------------------
%           Plate
% -------------------------
%                       C

% Material properties of the plate and dual PZTs size used in FEM
% simulation is identifical to the one used in numerical simulation for
% computing theoretical scaling matrix. "CompTheoreticalScalingMatirx"
% include all related properties

%% Load FEM simulation data
%
% We model measureant at a PZT as averaging of strains in X and Z
% directions (in-plane direction). We sum up those response in elements 
% corresponding to PZT areas (here, whole (ring+circle), ring, or circle).
% We do not include a code to convert FEM output to signals. Response
% signals are stored in a "FEM_DATA' file.
% 
% AB_F{i,j} indicates excitation at dual PZT A and sensing at dual PZT B. 
% i and j denote the different part(s) of the dual PZT activated for
% excitation and sensing. 1, 2, and 3 represent the entire, ring and
% circular parts of the dual PZT, respectively. For example, AB_F13 means
% the response signals measured at the circular part of the dual PZT B when
% both ring and circular part of dual PZT A is actuated. AC_F{i,j} is also
% defined as a similar manner but sensing at dual PZT C. 

% FEM_DATA contain AB_Fij and AC_Fij and a total of variables are 2 and
% size of each variable is a 3x3 cell array.

%% Parameters
clear; clc; close all; format shortg; warning off;

%% Setup folders
folderDataFEM      = fullfile(cd(cd('..')),'data','fem');

% Ouput folder
folderOut          = fullfile(cd(cd('..')),'output');

load(fullfile(folderDataFEM,'FEM_DATA.mat'));
%% Decomposition of S0 and A0 signals using collocated dual PZTs B and C
%
% In Section 3.2, S0 and A0 mode signals are obatined by adding and
% subtracting respones from dual PZTs B and C, respectively. (Please refer
% the following paper : S.B. Kim, H. Sohn, Instantaneous reference-free 
% crack detection based on polarization characteristics of piezoelectric 
% materials, Smart Mater. Struct. 16 (2007) 2375–2387.

% decomposition of S0 and A0 mode signals using collocate dual PZTs
% We will use cell arrays for simplification.
S0_F_Col = cellfun(@(x,y) plus(x,y)/2, AB_F, AC_F,'UniformOutput', 0);
A0_F_Col = cellfun(@(x,y) minus(x,y)/2, AB_F, AC_F,'UniformOutput', 0);

%% Computation of scaling factors from S0 and A0 signals (FEM)
%
% scaling factor obtained from FEM simulation 
[~, S_LOC] = cellfun(@(x) max(abs(x)),S0_F_Col);
[~, A_LOC] = cellfun(@(x) max(abs(x)),A0_F_Col);

scaleFact_S_FEM = cellfun(@(x,y) x(y),S0_F_Col,num2cell(S_LOC));
scaleFact_A_FEM = cellfun(@(x,y) x(y),A0_F_Col,num2cell(A_LOC));

% scaling matrix from FEM simulation
scalingMatrixFEM    = [reshape(scaleFact_S_FEM',9,1)/scaleFact_S_FEM(1) ...
                       reshape(scaleFact_A_FEM',9,1)/scaleFact_A_FEM(1)];
clearvars S_LOC A_LOC scaleFact_S_FEM scaleFact_A_FEM                   

%% Computation of theoretical scaling factors from S0 and A0 signals

% scaling matrix from numerical (theoretical) simulation
scalingMatrixTheory = CompTheoryScalingMatrix;

%% Comparison of scaling factors from FEM and numerical (theoretical) simulation
% 
y_S0(:,1) = scalingMatrixTheory(:,1);
y_S0(:,2) = scalingMatrixFEM(:,1); 

y_A0(:,1) = scalingMatrixTheory(:,2);
y_A0(:,2) = scalingMatrixFEM(:,2); 

figure(1)
h = stem(1:9,y_S0(1:9,:),'fill','--','linewidth',2);
legend('Theory', 'FEM');legend('location','northEast','orientation', 'horizontal');legend('boxoff');
set(h(1),'MarkerFaceColor','red','Marker','s','markerEdgecolor','none')
set(h(2),'MarkerFaceColor','blue','Marker','d','markeredgecolor','none')
set(gca, 'xTickLabel',{'S_1_1';'S_1_2';'S_1_3';'S_2_1'; ...
    'S_2_2';'S_2_3';'S_3_1';'S_3_2';'S_3_3';})
xlim([0 10]);ylim([0 2.5]);
set(gcf,'pos',[50 50 450 250]);
set(gca,'fontsize',12,'linewidth',2,'fontweight','bold','XTick' ,[1:9])
xlabel('\bf Normalized scaling factor');
ylabel('\bf Normalized amplitdue');
print(fullfile(folderOut,'FEM_ScalingS0'),'-djpeg','-r0')

figure(2)
h = stem(1:9,y_A0(1:9,:),'fill','--','linewidth',2);
legend('Theory', 'FEM');legend('location','northEast','orientation', 'horizontal');legend('boxoff');
set(h(1),'MarkerFaceColor','red','Marker','s','markeredgecolor','none')
set(h(2),'MarkerFaceColor','blue','Marker','d','markeredgecolor','none')
set(gca, 'xTickLabel',{'A_1_1';'A_1_2';'A_1_3';'A_2_1'; ...
    'A_2_2';'A_2_3';'A_3_1';'A_3_2';'A_3_3';})
xlim([0 10]);ylim([-8 8]);
set(gcf,'pos',[50 50 450 250]);
set(gca,'fontsize',12,'linewidth',2,'fontweight','bold','XTick' ,[1:9])
xlabel('\bf Normalized scaling factor');
ylabel('\bf Normalized amplitdue');
print(fullfile(folderOut,'FEM_ScalingA0'),'-djpeg','-r0')

%% Mode decomposition using dual PZTs
% S0 and A0 modes can be reconstred by mutiplying scaling factors and
% common factors Ms = S11*C(1) for S0 mode and Ma = A11*C(2) for A0 mode.
% We know V matrix (signal matrix) and scaling matrix (scaled with S11 and
% A11) so we can compute common factor Ms and Ma.

inv_S   = pinv(scalingMatrixTheory);

Ms      =  inv_S(1,:) * cell2mat(reshape(AB_F',9,1));
Ma      =  inv_S(2,:) * cell2mat(reshape(AB_F',9,1));


% FEM simulation parameter

ini_shift   =  250;         % intial shift due to toneburst
fs          =  5000000;     % sampling frequency

x_time = ((1:901) - ini_shift)/fs*1000;     % for ms

sT = scalingMatrixTheory;   % for simplicity

figure(3)
plot(x_time,sT(3,1)*Ms*10^6,'r',x_time, S0_F_Col{1,3}*10^6,':b','linewidth', 2);
% title('S0 modes in V13');
legend('\bf Proposed','\bf Collocated');
legend('location','northwest','orientation','horizontal');
legend('boxoff');
xlabel('\bf Time (ms)');ylabel('\bf Amplitude');
set(gcf,'pos',[50 50 450 250]);set(gca,'fontsize',10,'linewidth',2,'fontweight','bold')
xlim([0 1.3e-01]);
print(fullfile(folderOut,'FEM_S0_V13'),'-djpeg','-r0')

figure(4)
plot(x_time,sT(3,2)*Ma*10^6,'r',x_time, A0_F_Col{1,3}*10^6,':b','linewidth', 2);
% title('A0 modes in V13');
legend('\bf Proposed','\bf Collocated');
legend('location','northwest','orientation','horizontal');
legend('boxoff');
xlabel('\bf Time (ms)');ylabel('\bf Amplitude');
set(gcf,'pos',[50 50 450 250]);
set(gca,'fontsize',10,'linewidth',2,'fontweight','bold')
xlim([0 1.3e-01]);
print(fullfile(folderOut,'FEM_A0_V13'),'-djpeg','-r0')

figure(5)
plot(x_time,sT(8,1)*Ms*10^6,'r',x_time, S0_F_Col{3,2}*10^6,':b','linewidth', 2);
% title('S0 modes in V23');
legend('\bf Proposed','\bf Collocated');
legend('location','northwest','orientation','horizontal');
legend('boxoff');
xlabel('\bf Time (ms)');ylabel('\bf Amplitude');
set(gcf,'pos',[50 50 450 250]);
set(gca,'fontsize',10,'linewidth',2,'fontweight','bold')
xlim([0 1.3e-01]);
print(fullfile(folderOut,'FEM_S0_V32'),'-djpeg','-r0')

figure(6)
plot(x_time,sT(8,2)*Ma*10^6,'r',x_time, A0_F_Col{3,2}*10^6,':b','linewidth', 2);
% title('A0 modes in V23');
legend('\bf Proposed','\bf Collocated');
legend('location','northwest','orientation','horizontal');
legend('boxoff');
xlabel('\bf Time (ms)');ylabel('\bf Amplitude');
set(gcf,'pos',[50 50 450 250]);
set(gca,'fontsize',10,'linewidth',2,'fontweight','bold')
xlim([0 1.3e-01]);
print(fullfile(folderOut,'FEM_A0_V32'),'-djpeg','-r0')





