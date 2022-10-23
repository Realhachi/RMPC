%*************************************************
% Design of a state-feedback tracking controller
% for the two-tank system
%*************************************************

clear all;
clc;

%%
%*******************************************
% System matrices of the linearized model
% around:
% h1* = 0.8 m
% h2* = 0.4 m
% q1* = 0.0050596 m^3/s
% q2* = 0.0063246 m^3/s
%*******************************************

A = [-0.0891 0.0891
      0.0891 -0.2895];
  
B = [14.0845 0
     0 14.0845];

C=[1 0
   0 1];

D=[ 0 0
    0 0];

disp('Linear model of the plant:');
model = ss( A,B,C,D)

%************************
%Checking the stability
%************************
fprintf('\n\nPoles:\n');
disp(eig (A))

%****************************************************
% Checking if the system is controlable or not
%****************************************************
CO = ctrb(A,B);
disp('Rank of the controllability matrix:');
rank(CO)

%%
%******************************************************
% Computation of the feedback gain via pole-placement
%******************************************************

% Damping Ratio: 0.8
% Settling time: 20 seconds

dr = 0.8;
ts = 20;
wn = 4.6/(dr*ts);
alpha = -dr*wn;
beta = wn*sqrt(1-dr^2);

%desired closed loop poles
cl_poles=[alpha+beta*i alpha-beta*i]'

%Computing K

fprintf('\nThe state feedback gain is:\n');
K = place(A,B,cl_poles)
%Be aware: No eigenvalue should have a multiplicity greater than the 
%number of inputs.

%%
%*********************************************
% Reference Input - full state feedback
% Computation of the matrices Nx and Nu
%*********************************************
    
%number of states     
nx = 2;
%number of inputs
nu =2;
%number of outputs
ny=2;

big_A = [A B
         C D];

big_Y =[ zeros(nx,ny)
         eye(ny,ny) ];

big_N = inv(big_A)*big_Y;

fprintf('\nReference-Input  full state feedback. The matrices Nx and Nu are:');

Nx = big_N(1:nx,:)
Nu = big_N (nx+1:end,:)

%return

%%
%**************************
%  Integral control
%  We add two integrators
%**************************

%Constructing the Augmented system
NA = [ zeros(ny,ny)  C
       zeros(nx,ny)  A];

NB = [ D 
       B];
  

%checking the controlabillity of the Augmented system
disp('Rank of the controllability matrix of the augmented system:');
rank(ctrb(NA,NB))

%computing the feedback matrix of the augmented system
full_K = place(NA,NB,[cl_poles;-2; -2.1] );

fprintf('\nThe state feedback gains of the augmented system are:\n');

Ki = full_K(:,1:2)
Ks = full_K(:,3:end)