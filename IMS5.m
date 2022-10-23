%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author:Hachi Wu
%Date: 18/04/2022
%
%This is simulation of a numerical example from the Paper Robust delay dependent iterative learning
%fault-tolerant control for batch processes with state delay and actuator failures

%Installation package to be installed -- Yalmip, mosek
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
yalmip('clear');
close all
clc

del_umax = 1; % maximum input
del_umin = -1; % minimum input
del_ymin = -2.3;
del_ymax = 2.3;
dm = 1;
dM = 3;
%Model parameters
% Ac = [1.582 -0.5916;1 0];
% A_d = [-0.21 0.148; 0 0];
% Bc = [1;0];
% Cc = [1.69 1.419];
Ac = [0.985 0.0107;0.0078 0.9784];
A_d = [0.1057 0.0004; 0.0002 0.0807];
Bc = [64.4453;0.2559];
Cc = [1 0];
dt = 0.1;
[A,B,C,~] = c2dm(Ac,Bc,Cc,0,dt);

%Model data
N = diag([0.1 0.1]);
H = diag([0.1 0.2]);
H_d = diag([0.1 0.3]);
H_b = [0.2;0.1];
N_hat = [N; C*N];
H_hat = [H zeros(2,1)];
H_dhat = [H_d zeros(2,1)];
H_bhat = H_b;
E_hat = [0 0 1];
nx = 2;     % Number of states
nu = 1;     % Number of inputs
% MPC data
% Q1   = diag([5 2 1]);
Q1 = diag([1 1 1]);
R1 = 1;
delay_bound = 3;
% T = 199+delay_bound;% iteration times
T = 29+delay_bound;
% 定义变量
delta_u = sdpvar(repmat(nu,1,T),repmat(1,1,T));
delta_x = sdpvar(repmat(nx,1,T+1),repmat(1,1,T+1));
delta_y = sdpvar(repmat(nu,1,T+1),repmat(1,1,T+1));
y = sdpvar(repmat(nu,1,T+1),repmat(1,1,T+1));
c = sdpvar(repmat(nu,1,T+1),repmat(1,1,T+1));   %setpoint
e = sdpvar(repmat(nu,1,T+1),repmat(1,1,T+1));
% x_hat = [delta_x;e];          %也许需要预先申请空间
%Initial state
for i = 1:1:delay_bound
    delta_x{i} = [0.3;-0.5];
    delta_u{i} = 0.2; 
    delta_y{i} = 0;
    c{i} = 0.1;
    y{i} = 0;
    e{i} = y{i} - c{i};
    x_hat{i} = [delta_x{i}; e{i}];
%     y_d{i} = 15;
%     y_r{i} = y_d{i} + y{i};
%     y_r{i} = 0;
end

%controller
LMI1 = cell(16,16);
LMI2 = cell(4,4);
LMI3 = cell(4,4);
LMI4 = cell(4,4);
% LMI5 = cell(3,3);
% LMI6 = cell(3,3);

theta= sdpvar(1);
P1 = sdpvar(3,3,'symmetric');
T1 = sdpvar(3,3,'symmetric');
% M1 = sdpvar(3,3,'symmetric'); 
G1 = sdpvar(3,3,'symmetric');
L1 = sdpvar(3,3,'symmetric');
S1 = sdpvar(3,3,'symmetric');
S2 = sdpvar(3,3,'symmetric');
% M3 = sdpvar(3,3,'symmetric');
% M4 = sdpvar(3,3,'symmetric');
X1 = sdpvar(3,3,'symmetric');
X2 = sdpvar(3,3,'symmetric');
X3 = sdpvar(3,3,'symmetric');
Y1 = sdpvar(1,3);
D = dM*eye(3);
%D2 = dM*eye(3);
% 中间变量定义
Psi = P1+dM*T1+dM^2*((1+dM)/2)*G1;
phi = -L1+S2-X2;

% f = theta>=0;
% f = f + [LMI2>=0] + [LMI3>=0] + [LMI4>=0];
% f = f + [LMI2>=0] + [LMI3>=0] + [LMI5>=0] +[LMI6>=0];
%data save
delu_record = [];
delx_record = [];
y_record = [];
c_record = [];
e_record = [];
%x_d_record = [];
d_record = [];  %时滞
delu_record = [delu_record value(delta_u{delay_bound})]; 
delx_record = [delx_record value(delta_x{delay_bound})];
y_record = [y_record value(delta_y{delay_bound})];

for k = delay_bound+1:T
    d = unidrnd(delay_bound);
    d_record = [d_record d];
    F1 = cos(10*k);
    F2 = (sin(10*k)+1)/2;
    % uncertain perturbations
    delta_a = N*diag([F1 F2])*H;  
    delta_d = N*diag([F1 F2])*H_d;
    delta_b = N*diag([F1 F2])*H_b;
    delta_ah = N_hat*diag([F1 F2])*H_hat; 
    delta_dh = N_hat*diag([F1 F2])*H_dhat; 
    delta_bh = N_hat*diag([F1 F2])*H_bhat; 
    A1 = A + delta_a;
    A1_d = A_d +delta_d;
    B1 = B + delta_b;
    A_hat = [A zeros(2,1);C*A eye(1)];
    A_hat = A_hat + delta_ah;
    A_dhat = [A_d zeros(2,1);C*A_d eye(1)];
    A_dhat = A_dhat + delta_dh;
    B_hat = [B; C*B];
    B_hat = B_hat + delta_bh;
    C_hat = [C zeros(1,1)];
%     Phi = theta*inv(Psi);       %有问题，没有定义对Psi的逆
    x_hat{k} = [delta_x{k};e{k}];
    
    LMI1 = [phi                  L1          L1*A_hat'+Y1'*B_hat'       L1*A_hat'+Y1'*B_hat'-L1 Y1'*sqrt(R1) L1*sqrt(Q1);
            L1                   -X3-X1      X1*A_dhat'                 X1*A_dhat'              zeros(3,1)   zeros(3,3);
            A_hat*L1+B_hat*Y1    A_dhat*X1   -L1                        zeros(3,3)              zeros(3,1)   zeros(3,3);
            A_hat*L1+B_hat*Y1-L1 A_dhat*X1   zeros(3,3)                 -D^(-2)*X1              zeros(3,1)   zeros(3,3);        
            sqrt(R1')*Y1         zeros(1,3)  zeros(1,3)                 zeros(1,3)              -theta*eye(1) zeros(1,3);
            sqrt(Q1)*L1          zeros(3,3)  zeros(3,3)                 zeros(3,3)              zeros(3,1)   -theta*eye(3) ];    %稳定性约束1

%     LMI2 = [-1 x_hat{k}';
%             x_hat{k} -Phi];
    LMI2 = [-1  x_hat{k}'*Psi;
            Psi*x_hat{k} -theta*Psi];
%     LMI3 = [-del_umax^2 Y1
%             Y1'         -Phi];
    LMI3 = [-del_umax^2 Y1*Psi;
            Psi*Y1'     -theta*Psi];
%     LMI4 = [-del_ymax^2* C_hat;
%             C_hat'               -eye(1)];  %有争议，Phi的维数是几？？
    LMI4 = [-del_ymax^2*theta*Psi theta*C_hat';
            C_hat*theta           -1];
    f = [LMI1<=0] + [LMI2<=0] + [LMI3<=0] + [LMI4<=0];
%    f = [LMI1<=0] + [LMI2<=0] + [LMI3<=0];
    f = [f,del_umin <= delta_u{k} <= del_umax, del_ymin <= delta_y{k} <= del_ymax];
    
    ops = sdpsettings('warning',1,'verbose',1,'solver','mosek','cachesolvers',1);
    obj = theta;
    sol = optimize(f,obj);
    if sol.problem == 0
        disp('Solver thinks it is feasible');
    end
    Y_ = value(Y1);
    L_ = value(L1);
%     e{k} = y{k} - c{k};
    
    K_hat = Y_*inv(L_);
    delta_u{k} = K_hat*x_hat{k};
    
%     del_x{k+1} = A1{k}*del_x{k} + A1_d{k}*del_x{k-d} +B1{k}*delta_u{k};
%     delta_y{k} = C*delta_x{k};
%     x_hat{k+1} = A_hat{k}*x_hat{k} + A_dhat*x_hat{k-d} +B_hat{k}*delta_u{k};
%     x_hat{k+1} = (A_hat+B_hat*K_hat)*x_hat{k} + A_dhat*x_hat{k-d};
    x_hat{k+1} = (A_hat+B_hat*K_hat)*x_hat{k}+A_dhat*x_hat{k-d};
    delta_y{k} = C_hat*x_hat{k};
    e{k} = E_hat*x_hat{k};
    
    if k <= 42
        c{k} = 0.1;
%         y_d{k+1} = 15;
%         y_r{k+1} = y_d{k+1} + del_y{k+1};
    else
        c{k} = 0.2;
%         y_d{k+1} = 30;
%         y_r{k+1} = y_d{k+1} + del_y{k+1};
    end
    delu_record = [delu_record value(delta_u{k})];
    delx_record = [delx_record value(x_hat{k+1}(1:2))];
    y_record = [y_record value(delta_y{k})];
end    
  
figure
plot(delx_record(1,:),'linewidth',1);
hold on
plot(delx_record(2,:),'r--','linewidth',1);
xlabel('t/sampling')
ylabel('x(t)')
legend('x1','x2');

figure
plot(delu_record,'linewidth',1);
xlabel('t/sampling')
ylabel('u(t)')   

figure
plot(y_record,'r--','linewidth',1);
hold on 
plot(y_record,'linewidth',1);
xlabel('step')
ylabel('Output-Y')
legend('setpoint','Robust output');
