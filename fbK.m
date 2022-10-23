%%---------Record---------
%完成了反馈加跟踪控制
%%
clearvars
yalmip('clear');
close all
clc
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
dt = 0.1;

K = sdpvar(1,3);
P = sdpvar(3,3,'symmetric');
Q = sdpvar(3,3,'symmetric');
S = sdpvar(1,3);
% alfa = 10;

L1 = [-eye(3) Q*Az'+S'*Bz';Az*Q+Bz*S -Q]; 

L = [L1<=0, P>=0, Q>=0];
solvesdp(L);

% optimize(f,[],sdpsettings('solver','mosek'));
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
c = ones(1,T+1);   %setpoint
x = cell(1,T+1);
u = cell(1,T+1);
x{1} = [1;1];
u{1} = 0;
y{1} = 0;

delta_x{1} = [0.1;0.1];
delta_u{1} = 0.5; 
delta_y{1} = 0.5;
c(1) = 10;
%data save
delu_record = [];
z_record = [];
dely_record = [];
c_record = [];
e_record = [];
d_record = [];  %Time Delays
x_record = [x{1}];
u_record = [u{1}];
y_record = [y{1}];
e{1} = y{1} - c(1);
z{1} = [delta_x{1}; e{1}];
z_record = [z{1}];
for k = 1:1:T
    c(k+1)=10;
    delta_u{k+1} = K*z{k};
    u{k+1} = u{k} + delta_u{k+1}; 
%     x{k+1} = x{k} + delta_x{k};
    
    z{k+1} = Az*z{k} + Bz*delta_u{k};
    x{k+1} = z{k}(1:2) + x{k};
    delta_y{k+1} = Cz*z{k+1};
    y{k+1} = y{k} + delta_y{k+1};
    z_record = [z_record z{k+1}];
    delu_record = [delu_record value(delta_u{k})];
    z_record = [z_record value(z{k})];
    dely_record = [dely_record value(delta_y{k+1})];
    x_record = [x_record x{k+1} ];
    u_record = [u_record u{k+1}];
    y_record = [y_record y{k+1}];
    e_record = [e_record value(e{k})];
end
%%
if 0
figure(1)
plot(z_record(1,:),'linewidth',1);
hold on
plot(z_record(2,:),'r--','linewidth',1);
xlabel('t/sampling')
ylabel('x(t)')
legend('x1','x2');

figure(2)
plot(delu_record,'linewidth',1);
xlabel('t/sampling')
ylabel('u(t)')   

figure(3)
plot(c,'linewidth',1);
hold on 
plot(dely_record,'r--','linewidth',1);
xlabel('step')
ylabel('Output-Y')
legend('setpoint','statefeedback');
end
if 1
figure
plot(x_record(1,:),'linewidth',1);
% plot(z_record(1,:),'linewidth',1);
hold on
plot(x_record(2,:),'r--','linewidth',1);
% plot(z_record(2,:),'linewidth',1);
xlabel('t/sampling')
ylabel('x(t)')
legend('x1','x2');

figure
% plot(delu_record,'linewidth',1);
plot(u_record,'linewidth',1)
xlabel('t/sampling')
ylabel('u(t)')   

figure
% plot(dely_record,'r--','linewidth',1);
plot(c,'linewidth',1);
hold on 
plot(y_record,'r--','linewidth',1);
xlabel('step')
ylabel('Output-Y')
legend('setpoint','statefeedback');
end
