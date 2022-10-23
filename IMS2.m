%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author:Hachi Wu
%Date: 18/04/2022
%
%This is simulation of a numerical example from the Paper Robust delay dependent iterative learning
%fault-tolerant control for batch processes with state delay and actuator failures

%Installation package to be installed -- Yalmip, penlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
yalmip('clear');
close all
clc

umax = 5; % maximum input
umin = -5; % minimum input
%Model parameters
Ac = [1.582 -0.5916;1 0];
Ad = [-0.21 0.148; 0 0];
Bc = [1;0];
Cc = [1.69 1.419];

tao = 0.2;
%Model data
H1 = [0.03;0.02];
H2 = H1;
E1 = [0 1];
E2 = E1;
E3 = 0.1;
nx = 2;     % Number of states
nu = 1;     % Number of inputs
% MPC data
Q1   = eye(2);
R   = 1;
delay_bound = 3;
N = 29+delay_bound;% iteration times

dt = 0.1;
[A,B,C,~] = c2dm(Ac,Bc,Cc,0,dt);
% A = [-4 1;-2 -5];
% B = [0;1];

u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
%x_d = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
%x_r = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
y = sdpvar(repmat(nu,1,N+1),repmat(1,1,N+1));
y_d = sdpvar(repmat(nu,1,N+1),repmat(1,1,N+1));
y_r = sdpvar(repmat(nu,1,N+1),repmat(1,1,N+1));

%Initial state
for i = 1:1:delay_bound
    x{i} = [0.5;-0.5];
    u{i} = 0;
    y{i} = 0;
%    x_d{i} = [0;0];
%    x_r{i} = x_d{i} + x{i};
    y_d{i} = 15;
%     y_r{i} = y_d{i} + y{i};
    y_r{i} = 0;
end

%controller
LMI1 = cell(5,5);
LMI2 = cell(5,5);
LMI3 = cell(10,10);
LMI4 = cell(3,3);
LMI5 = cell(3,3);
LMI6 = cell(3,3);

gamma = sdpvar(1);
Q = sdpvar(2,2,'symmetric');
W = sdpvar(2,2,'symmetric');
Y = sdpvar(1,2);
Lambda1 = sdpvar(1);
Lambda2 = sdpvar(1); 
% LMI2 = [0.5*inv(W)   E2'         Ad';
%         E2      Lambda2     zeros(1,2);
%         Ad      zeros(2,1)  Q-Lambda2*(H2*H2')];
LMI2 = [0.5*W   W*E2'           W*Ad';
        E2*W    Lambda2*eye(1)  zeros(1,2);
        Ad*W    zeros(2,1)      Q-Lambda2*(H2*H2')];    %稳定性约束1
LMI3 = [0.5*Q       (Q1^(0.5)*Q)'   (R^(0.5)*Y)'    Q           (E1*Q+E3*Y)'    (A*Q+B*Y)';
        Q1^(0.5)*Q  2*gamma*eye(2)  zeros(2,1)      zeros(2,2)  zeros(2,1)      zeros(2,2);
        R^(0.5)*Y   zeros(1,2)      2*gamma*eye(1)  zeros(1,2)  zeros(1,1)      zeros(1,2);
        Q           zeros(2,2)      zeros(2,1)      2*W         zeros(2,1)      zeros(2,2);
        E1*Q+E3*Y   zeros(1,2)      zeros(1,1)      zeros(1,2)  Lambda1*eye(1)  zeros(1,2);
        A*Q+B*Y     zeros(2,2)      zeros(2,1)      zeros(2,2)  zeros(2,1)      Q-Lambda1*(H1*H1')];
LMI4 = [umax^2*eye(1)   Y;
        Y'              Q];             %输入约束
LMI5 = [Q           H1*Lambda1;
        Lambda1*H1' Lambda1*eye(1)];    %稳定性约束2
LMI6 = [Q           H2*Lambda2;
        Lambda2*H2' Lambda2*eye(1)];    %稳定性约束3
f = [Q>=0,gamma>=0,W>=0,Lambda1>=0,Lambda2>=0];
f = f + [LMI2>=0] + [LMI3>=0] + [LMI4>=0] + [LMI5>=0] +[LMI6>=0];
% f = f + [LMI2>=0] + [LMI3>=0] + [LMI5>=0] +[LMI6>=0];
%data save
u_record = [];
x_record = [];
y_d_record = [];
y_r_record = [];
%x_r_record = [];
%x_d_record = [];
d_record = [];  %时滞
x_record = [x_record value(x{delay_bound})];
%x_d_record = [x_d_record value(x_d{delay_bound})];
%x_r_record = [x_r_record value(x{delay_bound})+value(x_d{delay_bound})];
y_d_record = [y_d_record value(y_d{delay_bound})];
y_r_record = [y_r_record value(y{delay_bound})];
for k = delay_bound+1:N
    d = unidrnd(delay_bound);
    d_record = [d_record d];
    F1 = cos(10*k);
    F2 = (sin(10*k)+1)/2;
    A = A + H1*F1*E1;
    B = B + H1*F1*E3;
    Ad = Ad + H2*F2*E2;
    xi3 = [x{k}',x{k-1}',x{k-2}',x{k-3}']';
    dia = [Q            zeros(2,2) zeros(2,2) zeros(2,2);
           zeros(2,2)   W          zeros(2,2) zeros(2,2)
           zeros(2,2)   zeros(2,2) W          zeros(2,2)
           zeros(2,2)   zeros(2,2) zeros(2,2) W]; 
    LMI1 = [1   xi3';   %保证无限时域鲁棒模型性能指标上界的最小化
            xi3 dia];
    f = f + [LMI1>=0];
    f = [f,umin <= u{k} <= umax];
    
    ops = sdpsettings('warning',1,'verbose',1,'solver','mosek','cachesolvers',1);
    obj = gamma;
    sol = optimize(f,obj);
    if sol.problem == 0
        disp('Solver thinks it is feasible');
    end
    Y_ = value(Y);
    Q_ = value(Q);
    u{k} = Y_*inv(Q_)*x{k};
    x{k+1} = A*x{k} + Ad*x{k-d} + B*u{k};
    y{k+1} = C*x{k};
    y_d{k+1} = 15;
    y_r{k+1} = y_d{k+1} + y{k+1};
    
    u_record = [u_record value(u{k})];
    x_record = [x_record value(x{k+1})];
    y_d_record = [y_d_record value(y_d{k+1})];
    y_r_record = [y_r_record value(y_r{k+1})];
end
figure
plot(x_record(1,:),'linewidth',1);
hold on
plot(x_record(2,:),'r--','linewidth',1);
xlabel('t/sampling')
ylabel('x(t)')
legend('x1','x2');

figure
plot(u_record,'linewidth',1);
xlabel('t/sampling')
ylabel('u(t)')   

figure
plot(y_d_record,'r--','linewidth',1);
hold on 
plot(y_r_record,'linewidth',1);
xlabel('step')
ylabel('Output-Y')
legend('setpoint','Robust output');
