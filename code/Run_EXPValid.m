%% Validation of the proposed technique using experimental data
%
%% Contact
%  Name  : Chul Min Yeum
%  Email : chulminy@gmail.com
%  Please contact me if you have a question or find a bug

%% Reference
%
% Yeum, Chul Min, Hoon Sohn, and Jeong Beom Ihn. “Lamb Wave Mode
% Decomposition Using Concentric Ring and Circular Piezoelectric
% Transducers.” Wave Motion 48.4 (2011): 358–370.

%% Test configuration
%  Please refer to the author's original paper

%% Parameters
clear; clc; close all; format shortg; warning off;

%% Setup folders
folderDataExp       = fullfile(cd(cd('..')),'data','exp');

%% Load raw signals
% CH0 : Dual PZT A
% CH1 : Dual PZT B
% CH2 : Dual PZT D
% CH3 : Dual PZT C
%
% The channel order is that only dual PZTs A, B, and D are used
% for the proposed technique, which are installed on a single side. Dual
% PZT C (CH3) is installed to validate decomposition results.

% AB_F{i,j} indicates excitation at dual PZT A and sensing at dual PZT B.
% i and j denote the different part(s) of the dual PZT activated for
% excitation and sensing. 1, 2, and 3 represent the entire, ring and
% circular parts of the dual PZT, respectively. For example, AB_F13 means
% the response signals measured at the circular part of the dual PZT B when
% both ring and circular part of dual PZT A is actuated.

% Notations are defined in similar fashion in case of FEM simulation.
% The raw signals are stored in AB_F, AD_F, CB_F and CD_F.

AB_F = cell(3,3);AD_F = cell(3,3);
CB_F = cell(3,3);CD_F = cell(3,3);
for ii=1:3
    for jj=1:3
        AB = ['180k_CH0CH1_F' int2str(ii) int2str(jj)];
        AD = ['180k_CH0CH2_F' int2str(ii) int2str(jj)];
        CB = ['180k_CH3CH1_F' int2str(ii) int2str(jj)];
        CD = ['180k_CH3CH2_F' int2str(ii) int2str(jj)];
        load(fullfile(folderDataExp,AB));
        load(fullfile(folderDataExp,AD));
        load(fullfile(folderDataExp,CB));
        load(fullfile(folderDataExp,CD));
        AB_F{ii,jj} = eval(['transpose(X' AB ');']);
        AD_F{ii,jj} = eval(['transpose(X' AD ');']);
        CB_F{ii,jj} = eval(['transpose(X' CB ');']);
        CD_F{ii,jj} = eval(['transpose(X' CD ');']);
    end
end; clearvars -except AB_F AD_F CB_F CD_F;

%% Decomposition of S0 and A0 mode signals using collocate dual PZTs
% These are considered as reference (ground-truth).

% In Section 3.2, S0 and A0 mode signals are obatined by adding and
% subtracting respones from dual PZTs B and C, respectively. (Please refer
% the following paper : S.B. Kim, H. Sohn, Instantaneous reference-free
% crack detection based on polarization characteristics of piezoelectric
% materials, Smart Mater. Struct. 16 (2007) 2375–2387.

% S0_D_F_Col: S0 mode extracted from a signal measured between path AD
% using collocated PZTs
% A0_D_F_Col: A0 mode extracted from a signal measured between path AD
% using collocated PZTs
S0_D_F_Col = cellfun(@(x,y) plus(x,y)/2,  AD_F, CD_F,'UniformOutput',0);
A0_D_F_Col = cellfun(@(x,y) minus(x,y)/2, AD_F, CD_F,'UniformOutput',0);

% S0_B_F_Col: S0 mode extracted from a signal measured between path AB
% using collocated PZTs
% A0_B_F_Col: A0 mode extracted from a signal measured between path AB
% using collocated PZTs
S0_B_F_Col = cellfun(@(x,y) plus(x,y)/2,  AB_F, CB_F,'UniformOutput',0);
A0_B_F_Col = cellfun(@(x,y) minus(x,y)/2, AB_F, CB_F,'UniformOutput',0);

%% Decomposition of S0 and A0 mode at path AD using the proposed method
S0_range = 5000:5600;
A0_range = 5600:5900;

[~, S_LOC] = cellfun(@(x) max(x(S0_range)),AD_F);
[~, A_LOC] = cellfun(@(x) max(x(A0_range)),AD_F);

% scaling factor obtained from experiment signals
scaleFact_S_AD = cellfun(@(x,y) x(y),AD_F,num2cell(S_LOC+S0_range(1)-1));
scaleFact_A_AD = cellfun(@(x,y) x(y),AD_F,num2cell(A_LOC+A0_range(1)-1));

scalingMatrix_AD_MD  = [reshape(scaleFact_S_AD,9,1)/scaleFact_S_AD(1) ...
                        reshape(scaleFact_A_AD,9,1)/scaleFact_A_AD(1)];

% S0 and A0 modes can be reconstred by mutiplying scaling factors and
% common factors Ms = S11*C(1) for S0 mode and Ma = A11*C(2) for A0 mode.
% We know V matrix (signal matrix) and scaling matrix (scaled with S11 and
% A11) so we can compute common factor Ms and Ma. Please see eq.(7).
inv_S_AD    = pinv(scalingMatrix_AD_MD);
Ms_AD       = inv_S_AD(1,:) * cell2mat(reshape(AD_F,9,1));
Ma_AD       = inv_S_AD(2,:) * cell2mat(reshape(AD_F,9,1));

% S0 and A0 mode extracted using the proposed mode decomposition method
S0_D_F_MD = cellfun(@(x,y) x*y, repmat({Ms_AD},3,3), ...
    num2cell(reshape(scalingMatrix_AD_MD(:,1),3,3)), ...
    'UniformOutput',0);
A0_D_F_MD = cellfun(@(x,y) x*y, repmat({Ma_AD},3,3), ...
    num2cell(reshape(scalingMatrix_AD_MD(:,2),3,3)), ...
    'UniformOutput',0);

clearvars S_LOC A_LOC scaleFact_S_AD scaleFact_A_AD S0_range A0_range
clearvars inv_S_AD Ms_AD Ma_AD

%% Decomposition of S0 and A0 mode at path AD using the proposed method
inv_S_AB    = pinv(scalingMatrix_AD_MD);
Ms_AB       = inv_S_AB(1,:) * cell2mat(reshape(AB_F,9,1));
Ma_AB       = inv_S_AB(2,:) * cell2mat(reshape(AB_F,9,1));

S0_B_F_MD = cellfun(@(x,y) x*y, repmat({Ms_AB},3,3), ...
    num2cell(reshape(scalingMatrix_AD_MD(:,1),3,3)), ...
    'UniformOutput',0);
A0_B_F_MD = cellfun(@(x,y) x*y, repmat({Ma_AB},3,3), ...
    num2cell(reshape(scalingMatrix_AD_MD(:,2),3,3)), ...
    'UniformOutput',0);

clearvars S_LOC A_LOC scaleFact_S_AB scaleFact_A_AB S0_range A0_range
clearvars inv_S_AB Ms_AB Ma_AB

%% Output Plot

% Ouput folder
folderOut          = fullfile(cd(cd('..')),'output');

%% Comparison of the normalized scaling factors of the S0 and A0 modes
% In theory, the normalized scaling factors should be identical among 
% these two paths as long as the sizes of the used dual PZTs are identical.
% Please see Fig. 11.

% Scaling factor computation at AD using collocate PZTS
S0_range = 5000:5600;
A0_range = 5600:5900;

[~, S_LOC] = cellfun(@(x) max(x(S0_range)),S0_D_F_Col);
[~, A_LOC] = cellfun(@(x) max(x(A0_range)),A0_D_F_Col);

scaleFact_S_D = ...
    cellfun(@(x,y) x(y),S0_D_F_Col,num2cell(S_LOC+S0_range(1)-1));
scaleFact_A_D = ...
    cellfun(@(x,y) x(y),A0_D_F_Col,num2cell(A_LOC+A0_range(1)-1));

scalingMatrix_AD_Col  = [reshape(scaleFact_S_D',9,1)/scaleFact_S_D(1) ...
                         reshape(scaleFact_A_D',9,1)/scaleFact_A_D(1)];

clearvars S_LOC A_LOC scaleFact_S_D scaleFact_A_D S0_range A0_range

% Scaling factor computation at AB using collocate PZTS
S0_range = 5000:5400;
A0_range = 5000:5500;

[~, S_LOC] = cellfun(@(x) max(x(S0_range)),S0_B_F_Col);
[~, A_LOC] = cellfun(@(x) max(x(A0_range)),A0_B_F_Col);

scaleFact_S_B = ...
    cellfun(@(x,y) x(y),S0_B_F_Col,num2cell(S_LOC+S0_range(1)-1));
scaleFact_A_B = ...
    cellfun(@(x,y) x(y),A0_B_F_Col,num2cell(A_LOC+A0_range(1)-1));

scalingMatrix_AB_Col   = [reshape(scaleFact_S_B',9,1)/scaleFact_S_B(1) ...
                          reshape(scaleFact_A_B',9,1)/scaleFact_A_B(1)];

clearvars S_LOC A_LOC scaleFact_S_B scaleFact_A_B S0_range A0_range

y_S0(:,1) = scalingMatrix_AD_Col(:,1);
y_S0(:,2) = scalingMatrix_AB_Col(:,1);

y_A0(:,1) = scalingMatrix_AD_Col(:,2);
y_A0(:,2) = scalingMatrix_AB_Col(:,2);

figure(1) ; x = 1:9;
h = stem(x,y_S0,'fill','--','linewidth',2);
legend('AB path', 'AD path');legend('location','northEast', ...
    'orientation', 'horizontal');legend('boxoff');
set(h(1),'MarkerFaceColor','red','Marker','s','markeredgecolor','none')
set(h(2),'MarkerFaceColor','blue','Marker','d','markeredgecolor','none')
set(gca, 'xTickLabel',{'S_1_1';'S_1_2';'S_1_3';'S_2_1'; ...
    'S_2_2';'S_2_3';'S_3_1';'S_3_2';'S_3_3';})
xlim([0 10]);%ylim([-8 8]);
set(gca,'fontsize',12,'linewidth',2,'fontweight','bold','XTick' ,[1:9]);
set(gcf,'pos',[50 50 450 250]);
xlabel('\bf Normalized scaling factor');
ylabel('\bf Normalized amplitdue');
print(fullfile(folderOut,'EXP_ScalingS0'),'-djpeg','-r0')

figure(2) ; x = 1:9;
h = stem(x,y_A0,'fill','--','linewidth',2);
legend('AB path', 'AD path');legend('location','northEast', ...
    'orientation', 'horizontal');legend('boxoff');
set(h(1),'MarkerFaceColor','red','Marker','s','markeredgecolor','none')
set(h(2),'MarkerFaceColor','blue','Marker','d','markeredgecolor','none')
set(gca, 'xTickLabel',{'A_1_1';'A_1_2';'A_1_3';'A_2_1'; ...
    'A_2_2';'A_2_3';'A_3_1';'A_3_2';'A_3_3';})
xlim([0 10]);%ylim([-8 8]);
set(gca,'fontsize',12,'linewidth',2,'fontweight','bold','XTick' ,[1:9]);
set(gcf,'pos',[50 50 450 250]);
xlabel('\bf Normalized scaling factor');
ylabel('\bf Normalized amplitdue');
print(fullfile(folderOut,'EXP_ScalingA0'),'-djpeg','-r0')


%% Comparison between the S0 and A0 modes in the path AD 
% decomposed by the proposed technique and the ones selectively generated
% by the collocated PTS. Please see Fig. 10.

% experiment parameter
startTime   =  5000;        % time point at a peak of toneburst input
fs          =  5000000;     % sampling frequency
endTime     =  6500;        % end time point to be plotted

x_path = 5000:6500;
x_time = (x_path-5000)./fs*1000 ;  % for ms

figure(3); 
plot(x_time, S0_B_F_Col{1,3}(x_path),'r' , x_time, ...
    S0_B_F_MD{1,3}(x_path),':b','linewidth',2)
set(gca,'Ytick',[-0.8 0 0.8]);
legend('Collocated','Proposed','Orientation','horizontal', ...
    'location','north');legend('boxoff')
set(gca,'Xtick',([5000 5300 5600 5900 6200 6500]-5000)./5000000*1000);
set(gca,'fontsize',12,'linewidth',1,'fontweight','bold')
xlim(([5000 6500]-5000)./5000000*1000); ylim([-0.8 0.8]);
ylabel('\bf Voltage(V)');xlabel('\bf Time(ms)');
set(gcf,'pos',[50 50 450 250]);
print(fullfile(folderOut,'EXP_AB_S0_V13'),'-djpeg','-r0');

figure(4);
plot(x_time, A0_B_F_Col{1,3}(x_path),'r' , x_time, ...
    A0_B_F_MD{1,3}(x_path),':b','linewidth',2)
set(gca,'Ytick',[-0.8 0 0.8]);
legend('Collocated','Proposed','Orientation','horizontal', ...
    'location','north');legend('boxoff')
set(gca,'Xtick',([5000 5300 5600 5900 6200 6500]-5000)./5000000*1000);
set(gca,'fontsize',12,'linewidth',1,'fontweight','bold')
xlim(([5000 6500]-5000)./5000000*1000); ylim([-0.8 0.8]);
ylabel('\bf Voltage(V)');xlabel('\bf Time(ms)');
set(gcf,'pos',[50 50 450 250]);
print(fullfile(folderOut,'EXP_AB_A0_V13'),'-djpeg','-r0');

figure(5);
plot(x_time, S0_B_F_Col{3,2}(x_path),'r' , x_time, ...
    S0_B_F_MD{3,2}(x_path),':b','linewidth',2)
set(gca,'Ytick',[-0.3 0 0.3]);
legend('Collocated','Proposed','Orientation','horizontal', ...
    'location','north');legend('boxoff')
set(gca,'Xtick',([5000 5300 5600 5900 6200 6500]-5000)./5000000*1000);
set(gca,'fontsize',12,'linewidth',1,'fontweight','bold')
xlim(([5000 6500]-5000)./5000000*1000); ylim([-0.3 0.3]);
set(gcf,'pos',[50 50 450 250]);
print(fullfile(folderOut,'EXP_AB_S0_V32'),'-djpeg','-r0');

figure(6);
plot(x_time, A0_B_F_Col{3,2}(x_path),'r' , x_time, ...
    A0_B_F_MD{3,2}(x_path),':b','linewidth',2)
set(gca,'Ytick',[-0.3 0 0.3]);
legend('Collocated','Proposed','Orientation','horizontal', ...
    'location','north');legend('boxoff')
set(gca,'Xtick',([5000 5300 5600 5900 6200 6500]-5000)./5000000*1000);
set(gca,'fontsize',12,'linewidth',1,'fontweight','bold')
xlim(([5000 6500]-5000)./5000000*1000); ylim([-0.3 0.3]);
ylabel('\bf Voltage(V)');xlabel('\bf Time(ms)');
set(gcf,'pos',[50 50 450 250]);
print(fullfile(folderOut,'EXP_AB_A0_V32'),'-djpeg','-r0');