%%---------Record---------
%目标完成MPC求反馈增益K
%%
clearvars
yalmip('clear');
close all
clc
del_umax = 10; % maximum input
delta_ymax = 15;
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
Gp = [eye(2);C];
E = [0.02;0.4];     %disturbance related coeffi.


% MPC data
Q1 = diag([30,15,20]);
R = 1;
rou = 3;
beta = 1;   % fault related coeffi.
nx = 2;     % Number of states
nu = 1;     % Number of inputs
delay_bound = 2;
T = 59+delay_bound;
delta_u = cell(1,T);
delta_x = cell(2,T+1);
% delta_x = sdpvar(repmat(nx,1,T+1),repmat(1,1,T+1));
delta_y = cell(1,T+1);
y = cell(1,T+1);
e = cell(1,T+1);
% e = sdpvar(repmat(nx,1,T+1),repmat(1,1,T+1));
c = cell(1,T+1);   %setpoint
x = cell(1,T+1);
u = cell(1,T+1);
x{1} = [1;1];
u{1} = 0;
y{1} = 0;
delta_x{1} = [0.1;0.1];
delta_u{1} = 0.5; 
delta_y{1} = 0.5;
c{1} = 10;
%data save
delu_rec = [];
z_rec = [];
dely_rec = [];
c_rec = [];
% e_rec = [];
d_rec = [];  %Time Delays
x_rec = [x{1}];
u_rec = [u{1}];
y_rec = [y{1}];
F_rec = [];

e{1} = y{1} - c{1};
z{1} = [delta_x{1}; e{1}];
z_rec = [z{1}];

Q = sdpvar(3,3,'symmetric');
% W = sdpvar(3,3,'symmetric');
Y = sdpvar(1,3);
gamma = sdpvar(1);
%Initial state

%controller
% LMI1 = cell(16,16);
% LMI2 = cell(4,4);
% LMI3 = cell(4,4);
% LMI4 = cell(4,4);
LMI2 = [Q           Q*Az'+Y'*Bz' Q*sqrt(Q1)     Y'*sqrt(R);
        Az*Q+Bz*Y   Q            zeros(3,3)     zeros(3,1);
        sqrt(Q1)*Q  zeros(3,3)   gamma*eye(3)   zeros(3,1);
        sqrt(R)*Y   zeros(1,3)   zeros(1,3)     gamma*eye(1)];
LMI3 = [del_umax^2*eye(1)   Y;
        Y'                  Q];    %Input constraint
LMI4 = [Q               (Az*Q+Bz*Y)'*Cz';
        Cz*(Az*Q+Bz*Y)  delta_ymax^2*eye(1)];   %output constraint
f = [Q>=0,gamma>=0];
f = f + [LMI2>=0];
for k = 2:1:T+1
    c{k}=10;
    LMI1 = [1       z{k-1}';
            z{k-1}    Q   ];   %/tao = 2 
    f = f + [LMI1>=0];
   
    ops = sdpsettings('warning',1,'verbose',1,'solver','sedumi','cachesolvers',1);
    obj = gamma;
    sol = optimize(f,obj);
    if sol.problem == 0
        disp('Solver thinks it is feasible');
    end
    Y_ = value(Y);
    Q_ = value(Q);
    F = Y_*inv(Q_);
    delta_u{k} = F*z{k-1};
    u{k} = u{k-1} + delta_u{k}; 
%     x{k+1} = x{k} + delta_x{k};
    z{k} = Az*z{k-1} + Bz*delta_u{k};
    x{k} = z{k}(1:2) + x{k-1};
    delta_y{k} = Cz*z{k};
    y{k} = y{k-1} + delta_y{k};
    z_rec = [z_rec z{k}];
    delu_rec = [delu_rec value(delta_u{k})];
    dely_rec = [dely_rec value(delta_y{k})];
    x_rec = [x_rec x{k} ];
    u_rec = [u_rec u{k}];
    y_rec = [y_rec y{k}];
%     e_rec = [e_rec value(e{k})];
    F_rec = [F_rec; F];
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
