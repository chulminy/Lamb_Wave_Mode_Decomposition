%% Theoretical computation of a scaling matrix using analytical solution of the Lamb wave response at a circular PZT
%
%%
function scalingMatrixTheory = CompTheoryScalingMatrix

%% Contact
%  Name  : Chul Min Yeum
%  Email : chulminy@gmail.com
%  Please contact me if you have a question or find a bug

%% Reference
% 
% * Yeum, Chul Min, Hoon Sohn, and Jeong Beom Ihn. “Lamb Wave Mode
% Decomposition Using Concentric Ring and Circular Piezoelectric
% Transducers.” Wave Motion 48.4 (2011): 358–370.
% * Giurgiutiu, Victor. Structural health monitoring: with piezoelectric 
% wafer active sensors. Academic Press, 2007.
% * Raghavan, Ajay, and Carlos ES Cesnik. "Modeling of piezoelectric-based
% Lamb wave generation and sensing for structural health monitoring." Smart
% Structures and Materials. International Society for Optics and Photonics,
% 2004.
% * Sohn, Hoon, and Sang Jun Lee. "Lamb wave tuning curve calibration for
% surface-bonded piezoelectric transducers." Smart Materials and Structures
% 19.1 (2009): 015007.
% * Giurgiutiu, Victor. Structural health monitoring: with piezoelectric
% wafer active sensors. Academic Press, (Chapter 6: Guided Waves),2007.

%% Variables
% Naming convention of variables in this code are close to notations used 
%
% * h           :   thickness of the plate structure (m)
% * b           :   h/2 (m) 
% * wb          :   width of the place (m)
% * E           :   Young's modulus of the plate (Pa)
% * nu          :   Poisson ratio of the plate
% * rho         :   density [kg/m^3]
% * Lamda       :   Lame' constant
% * Mu          :   Lame' constant
% * cp          :   P-wave velocity [m/s]
% * cs          :   S-wave velocity [m/s]
% * a1          :   radius of the outer ring actuator       [m] 
% * a2          :   radius of the inner ring actuator       [m]
% * a3          :   radius of the inner cicular actuator    [m]
% * c1          :   radius of the outer ring sensor         [m] 
% * c2          :   radius of the inner ring sensor         [m]
% * c3          :   radius of the inner cicular sensor      [m]
% * rs          :   distance between an actuator and the second sensor [m]
% * dr          :   spatial resolution [m]
% * xi_S0       :   wavenumber of S0 mode
% * xi_A0       :   wavenumber of A0 mode
% * freq_D      :   driving frequency

%% Initialize parameters

% aluminium plate(T-6061)
h       = 0.003;                % thickness of the plate [m]
E       = 7*1e10;               % Young modulus of the aluminum plate [Pa]
nu      = 0.3;                  % Poisson ratio of the aluminum plate
rho     = 2700;                 % density [kg/m^3]
Lamda   = E*nu/(1+nu)/(1-2*nu); % Lame' constant
Mu      = E/2/(1+nu);           % Lame' constant

cp      = sqrt((Lamda+2*Mu)/rho);    % P-wave velocity [m/s]
cs      = sqrt(Mu/rho);              % S-wave velocity [m/s]

% dimension of a dual PZT (actuar and sensor are same size)
a1      = 0.009;                % radius of the outer ring dual PZT     [m] 
a2      = 0.005;                % radius of the inner ring dual PZT     [m]
a3      = 0.004;                % radius of the inner cicular dual PZT  [m]

rs      = 0.2 ;                 % distance between an act. and the sens.[m]
dr      = 1e-6;                 % spatial resolution

freq_D  = 180;                  % driving frequency (kHz)

% calculate wavenumber
w = 2*pi*freq_D*1000;           % input frequency (kHz)

c0    = 1.8*cs;                      % intiatial wave speed (Figure 6.10)
xi_S0 = symme(c0,cp,cs,h,w);         % wavenumber for symmetric mode
% fd = freq_D*1000*h/2 = 270 (less than 1000). c/cs is around 1.8 in Fig.
% 6.10.

c0    = 0.5*cs;                      % intiatial wave speed (Figure 6.12)
xi_A0 = anti(c0,cp,cs,h,w);          % wavenumber for antisymmetric mode
% fd = freq_D*1000*h/2 = 270 (less than 1000). c/cs is around 1 in Fig.
% 6.12

%% Compute a scaling matrix 
%
K = pi * [a1*a1 a2*a2 a3*a3];   % Eq.(B4)

%%
% *S0 mode*

% temporaray variables
tA      = arrayfun(@(x) x*besselj(1,xi_S0*x),[a1 a2 a3]); 
tC      = arrayfun(@(x) besselj(1,xi_S0*x)/x,[a1 a2 a3]); 

S_S0 = repmat(tA',1,3).*repmat(tC,3,1); clearvars tA tC;   % Eq(6)
% Example: S_S0(1,2) indicates S^S0(a1,a2)  

% Eq.(B4)
S33 = S_S0(3,3);
S23 = S_S0(1,3)-S_S0(2,3);
S32 = (K(1)-K(2))^(-1)*(K(1)*S_S0(3,1)-K(2)*S_S0(3,2));
S22 = (K(1)-K(2))^(-1)* (...
      K(1)*(S_S0(1,1)-S_S0(2,1)) - K(2)*(S_S0(1,2)-S_S0(2,2)));
S21 = (K(1)-K(2)+K(3))^(-1) * ((K(1)-K(2))*S22 + K(3)*S23);
S31 = (K(1)-K(2)+K(3))^(-1) * ((K(1)-K(2))*S32 + K(3)*S33);

S13 = S23 + S33;
S12 = S22 + S32;
S11 = S21 + S31;

%%
% *A0 mode*

% temporaray variables
tA      = arrayfun(@(x) x*besselj(1,xi_A0*x),[a1 a2 a3]); 
tC      = arrayfun(@(x) besselj(1,xi_A0*x)/x,[a1 a2 a3]); 

S_A0    = repmat(tA',1,3).*repmat(tC,3,1); clearvars tA tC;   % Eq(6)

% Eq.(B4)
A33 = S_A0(3,3);
A23 = S_A0(1,3)-S_A0(2,3);
A32 = (K(1)-K(2))^(-1)*(K(1)*S_A0(3,1)-K(2)*S_A0(3,2));
A22 = (K(1)-K(2))^(-1)* (...
      K(1)*(S_A0(1,1)-S_A0(2,1)) - K(2)*(S_A0(1,2)-S_A0(2,2)));
A21 = (K(1)-K(2)+K(3))^(-1) * ((K(1)-K(2))*A22 + K(3)*A23);
A31 = (K(1)-K(2)+K(3))^(-1) * ((K(1)-K(2))*A32 + K(3)*A33);

A13 = A23 + A33;
A12 = A22 + A32;
A11 = A21 + A31;


%%
% *Scaling Matrix*
scalingMatrixTheory =  [  S11 A11; S12 A12; S13 A13; ...
                          S21 A21; S22 A22; S23 A23; ...
                          S31 A31; S32 A32; S33 A33];
scalingMatrixTheory(:,1) = scalingMatrixTheory(:,1)./S11;
scalingMatrixTheory(:,2) = scalingMatrixTheory(:,2)./A11;

end

%%
function xi = symme(c0,cp,cs,h,w)

% *Reference*
% * Giurgiutiu, Victor. Structural health monitoring: with piezoelectric
% wafer active sensors. Academic Press, (Chapter 6: Guided Waves), 2007.

% for symmetric case
% Solving the Rayleigh-Lamb frequcny relations to derive a wavenumber(xi)
options = optimset('Display','off');

d   = h/2;
f   = @(c) symme_fun(c,d,w,cp,cs);
xi  = w/fsolve(f,c0,options);

end

%%
function DS = symme_fun(c,d,w,cp,cs)
% *Reference*
% * Giurgiutiu, Victor. Structural health monitoring: with piezoelectric
% wafer active sensors. Academic Press, (Chapter 6: Guided Waves), 2007

xi = w/c;                                  % wave number
p  = sqrt((w/cp)^2 - xi^2);                % define 'p' coefficient Eq.82
q  = sqrt((w/cs)^2 - xi^2);                % define 'q' coefficient Eq.82
DS = (xi^2-q^2)^2*tan(q*d) + (4*xi^2*p*q)*tan(p*d); % Eq.101

end

%%
function xi = anti(c0,cp,cs,h,w)

% *Reference*
% * Giurgiutiu, Victor. Structural health monitoring: with piezoelectric
% wafer active sensors. Academic Press, (Chapter 6: Guided Waves), 2007.

% for anti symmetric case
% Solving the Rayleigh-Lamb frequcny relations to derive a wavenumber(xi)
options = optimset('Display','off');

d   = h/2;
f   = @(c) anti_fun(c,d,w,cp,cs);
xi  = w/fsolve(f,c0,options);

end

%%
function DA = anti_fun(c,d,w,cp,cs)

% *Reference*
% * Giurgiutiu, Victor. Structural health monitoring: with piezoelectric
% wafer active sensors. Academic Press, (Chapter 6: Guided Waves), 2007

xi = w/c;                                  % wave number
p  = sqrt((w/cp)^2 - xi^2);                % define 'p' coefficient Eq.82
q  = sqrt((w/cs)^2 - xi^2);                % define 'q' coefficient Eq.82
DA = (xi^2-q^2)^2*tan(p*d) + (4*xi^2*p*q)*tan(q*d); % Eq.110

end
