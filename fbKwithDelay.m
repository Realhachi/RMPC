%%---------Record---------
%目标完成带时滞的状态反馈跟踪控制
%%
clearvars
yalmip('clear');
close all
clc
% del_umax = 10; % maximum input
% del_umin = -10; % minimum input
A1 = [0.6703 0.0006;0.002 0.3679];
A2 = [0.5444 -0.0020;4.6488 -0.0736];
A3 = [0.4889 -0.0040;4.6488 -0.0736];
A4 = [0.4889 -0.0040;9.2976 0.0928];
B = [0.8242;0.0013];
C = [1 0];
N = diag([0.16 0.16]);
Az = [A2 zeros(2,1);C*A2 1];
Bz = [B;C*B];
Cz = [C 0];
Nz = [N zeros(2,1);C*N 0];
% Gp = [eye(2);C];
% E = [0.02;0.4];     %disturbance related coeffi.
% MPC data
% Q = diag([30,15,20]);
% R = 0.1;
% rou = 3;
beta = 1;   % fault related coeffi.

K = sdpvar(1,3);
P = sdpvar(3,3,'symmetric');
Q = sdpvar(3,3,'symmetric');
S = sdpvar(1,3);
L1 = [-eye(3) Q*Az'+S'*Bz';Az*Q+Bz*S -Q]; 
L = [L1<=0, P>=0, Q>=0];
solvesdp(L);
S_ = double(S);
Q_ = double(Q);
%statefeedback gain
K = S_*inv(Q_);

nx = 2;     % Number of states
nu = 1;     % Number of inputs
delay_bound = 2;
T = 59+delay_bound;
delta_u = cell(1,T);
delta_x = cell(2,T+1);
delta_y = cell(1,T+1);
y = cell(1,T+1);
e = cell(1,T+1);
c = cell(1,T+1);   %setpoint
x = cell(1,T+1);
u = cell(1,T+1);
x{1} = [1;1];
u{1} = 0;
y{1} = 0;

%data save
delu_rec = [];
z_rec = [];
dely_rec = [];
c_rec = [];
e_rec = [];
d_rec = [];  %Time Delays
x_rec = [x{1}];
u_rec = [u{1}];
y_rec = [y{1}];
%Initial state

for i = 1:1:delay_bound+1
    delta_x{i} = [0.1;0.1];
    delta_u{i} = 0.5;
    delta_y{i} = 0.5;
    x{i} = [1;1];
    u{i} = 0;
    y{i} = 1;
    c{i} = 10;
    e{i} = y{i} - c{i};
    x{i+1} = x{i} + delta_x{i};
    z{i} = [delta_x{i};e{i}];
    u{i+1} = u{i} + delta_u{i};
    y{i+1} = y{i} + delta_y{i};
    delu_rec=[delu_rec delta_u{i}];
    dely_rec = [dely_rec delta_y{i}];
    u_rec = [u_rec u{i+1}];
    x_rec = [x_rec x{i+1}];
    z_rec = [z_rec z{i}];
    e_rec = [e_rec e{i}];
    y_rec = [y_rec y{i+1}]; 
end
%controller
% LMI1 = cell(16,16);
% LMI2 = cell(4,4);
% LMI3 = cell(4,4);
% LMI4 = cell(4,4);
% LMI2 = [Omega11     zeros(3,3)   zeros(3,2)           G'          Omega51'    (sqrt(Q)*G)'    (sqrt(R)*Y)';
%         zeros(3,3)  W            zeros(3,2)           zeros(3,3)  (Nz*W)'     zeros(3,3)       zeros(3,1);
%         zeros(2,3)  zeros(2,3)   rou^2*gamma*eye(2)   zeros(2,3)  Gp'         zeros(2,3)       zeros(2,1);
%         G           zeros(3,3)   zeros(3,2)           W           zeros(3,3)  zeros(3,3)       zeros(3,1);
%         Omega51     Nz*W         Gp                   zeros(3,3)  Qlj         zeros(3,3)       zeros(3,1) ; 
%         sqrt(Q)*G   zeros(3,3)   zeros(3,2)           zeros(3,3)  zeros(3,3)  gamma*eye(3)     zeros(3,1);
%         sqrt(R)*Y   zeros(1,3)   zeros(1,2)           zeros(1,3)  zeros(1,3)  zeros(1,3)       gamma];
% LMI3 = [del_umax^2  Y;
%         Y'          Omega11];
for k = delay_bound+1:1:T
    d  = unidrnd(delay_bound); %产生随机数
    c{k+1}=10;
%     LMI1 = [1       z{k}'       z{k-1}'     z{k-2}';
%             z{k}    Qvj         zeros(3,3)  zeros(3,3);
%             z{k-1}  zeros(3,3)  W           zeros(3,3);
%             z{k-2}  zeros(3,3)  zeros(3,3)  W];   %/tao = 2 
%     f = [LMI1>=0, LMI2>=0, LMI3>=0];
%     f = [f, Qvj>=0, Qlj>=0, W>=0];
%     ops = sdpsettings('warning',1,'verbose',1,'solver','sedumi','cachesolvers',1);
%     obj = gamma;
%     sol = optimize(f,obj);
%     if sol.problem == 0
%         disp('Solver thinks it is feasible');
%     end
%     Y_ = value(Y);
%     G_ = value(G);
%     F = Y_*inv(G_);
    delta_u{k} = K*z{k};
    u{k+1} = u{k} + delta_u{k}; 
%     x{k+1} = x{k} + delta_x{k};
    z{k+1} = Az*z{k} + Bz*delta_u{k} + Nz*z{k-d};
    x{k+1} = z{k+1}(1:2) + x{k};
    delta_y{k+1} = Cz*z{k+1};
    y{k+1} = y{k} + delta_y{k+1};
    z_rec = [z_rec z{k+1}];
    delu_rec = [delu_rec value(delta_u{k})];
    dely_rec = [dely_rec value(delta_y{k+1})];
    x_rec = [x_rec x{k+1} ];
    u_rec = [u_rec u{k+1}];
    y_rec = [y_rec y{k+1}];
    e_rec = [e_rec value(e{k})];
end
%%
if 0
figure(1)
plot(z_rec(1,:),'linewidth',1);
hold on
plot(z_rec(2,:),'r--','linewidth',1);
xlabel('t/sampling')
ylabel('x(t)')
legend('x1','x2');

figure(2)
plot(delu_rec,'linewidth',1);
xlabel('t/sampling')
ylabel('u(t)')   

figure(3)
plot(c,'linewidth',1);
hold on 
plot(dely_rec,'r--','linewidth',1);
xlabel('step')
ylabel('Output-Y')
legend('setpoint','statefeedback');
end
if 1
figure
plot(x_rec(1,:),'linewidth',1);
% plot(z_rec(1,:),'linewidth',1);
hold on
plot(x_rec(2,:),'r--','linewidth',1);
% plot(z_rec(2,:),'linewidth',1);
xlabel('t/sampling')
ylabel('x(t)')
legend('x1','x2');

figure
% plot(delu_rec,'linewidth',1);
plot(u_rec,'linewidth',1)
xlabel('t/sampling')
ylabel('u(t)')   

figure
% plot(dely_rec,'r--','linewidth',1);
c = cell2mat(c);
plot(c,'linewidth',1);
hold on 
plot(y_rec,'r--','linewidth',1);
xlabel('step')
ylabel('Output-Y')
legend('setpoint','statefeedback');
end
