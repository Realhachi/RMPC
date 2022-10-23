%----With delays-----
clear all
%----- PARAMETERS DEFINITION -----
% A = [0.5444 -0.0020;2.3244 -0.1568];
% B = [0;0.8368];
% C = [0 1];
N = diag([0.16 0.16]);
% D = 0;
A  = [-3 1; 0 -2];
B = [0; 1];
C = [1 0];
D = [0];
%----- OPEN-LOOP ANALYSIS ----------------
openloop_eigenvalues=eig(A)
controlability_rank=rank(ctrb(A,B))
%----- STATE FEEDBACK DESIGN ------------
P=[-0.5+0.5i -0.5-0.5i];
% Desired Closed-Loop Ploes
K=acker(A,B,P); % State-Feedback Gain
%--------------------------------------------
%----- SIMULATION OF REGULATION PROBLEM -----
%--------------------------------------------
t=.01:.01:30; % simulation time
[rt,ct]=size(t);
for k = 1:1:4
    x(:,k)=[0;1]; % initial condition
end


% for k = 2:ct
%     xd(:,k-1) = [0;10];   %setpoint 
%     e(:,k-1) = xd(:,k-1)-x(:,k-1);
%     u(k) = -K*e(:,k-1);
% %     x(:,k) = x(:,k-1) + A*x(:,k-1)+B*u(k);
%     x(:,k) = x(:,k-1) + .01*(A*x(:,k-1)+B*u(k));
% end
%----------------------------------------------
%------ SIMULATION OF TRACKING PROBLEM --------
%----------------------------------------------
CL_sys=ss(A-B*K,B,C-D*K,D); % Closed-Loop System
CL_gain=dcgain(CL_sys);
% DC-Gain of closed-loop system
K_r=inv(CL_gain); % control signal: u=-K*x+K_r*yd
yd=10; % desired value for output (teta)
delay_bound = 2;
for k=4:ct
    d = unidrnd(delay_bound);
    u(k)=-K*x(:,k-1)+K_r*yd;
    x(:,k) = x(:,k-1)+.01*(A*x(:,k-1)+B*u(k)+N*x(:,k-d-1));
    y(:,k)=C*x(:,k-1);
end
figure;subplot(1,2,1);plot(t,x(1,:));
title('Regulation Problem, state:X_1');
subplot(1,2,2);plot(t,x(2,:));
title('Regulation Problem, state:X_2');
%~~~~~ PLOTTING CONTROL SIGNAL ~~~~
figure;plot(t,u);
title('Regulation Problem, Control Signal (u)');
%~~~~~ PLOTTING OUTPUT SIGNAL ~~~~
figure;plot(t,y);
title('Regulation Problem, Output Signal (y)');


