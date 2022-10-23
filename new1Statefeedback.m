clear all
%----- PARAMETERS DEFINITION -----
a11=-0.0507;
a12=-3.861;
a21=-0.00117;
a22=-0.5164;
a31=-0.000129;
a32=1.4168;
a33=-0.4932;
b1=0;
b2=-0.0717;
b3=-1.645;
g=9.81;
%----- OPEN-LOOP SYSTEM: X_dot=AX+BU ----
A=[a11 a12 0 -g;a21 a22 1 0;a31 a31 a33 0;0 1 0 0];
B=[b1;b2;b3;0];
%----- OPEN-LOOP ANALYSIS ----------------
openloop_eigenvalues=eig(A)
controlability_rank=rank(ctrb(A,B))
%----- STATE FEEDBACK DESIGN ------------
P=[-0.9+1.9i -0.9-1.9i -1.8+0.6i -1.8-.6i];
% Desired Closed-Loop Ploes
K=acker(A,B,P); % State-Feedback Gain
%--------------------------------------------
%----- SIMULATION OF REGULATION PROBLEM -----
%--------------------------------------------
t=.01:.01:10; % simulation time
[rt,ct]=size(t);
x(:,1)=[0;1;0;0]; % initial condition
for k=2:ct
u(k)=-K*x(:,k-1);
% x(:,k)=x(:,k-1)+.01*(A*x(:,k-1)+B*u(k));
x(:,k)=x(:,k-1)+1*(A*x(:,k-1)+B*u(k));
end
%~~~~~ PLOTTING STATES ~~~~~
subplot(2,2,1);plot(t,x(1,:));
title('Regulation Problem, state:X_1');
subplot(2,2,2);plot(t,x(2,:));
title('Regulation Problem, state:X_2');
subplot(2,2,3);plot(t,x(3,:));
title('Regulation Problem, state:X_3');
subplot(2,2,4);plot(t,x(4,:));
title('Regulation Problem, state:X_4');
%~~~~~ PLOTTING CONTROL SIGNAL ~~~~
figure;plot(t,u);
title('Regulation Problem, Control Signal (u)');
%*****************************
%----------------------------------------------
%------ SIMULATION OF TRACKING PROBLEM --------
%----------------------------------------------
C=[1 0 0 0]; % selecting 'teta' as system output
D=0;
CL_sys=ss(A-B*K,B,C-D*K,D); % Closed-Loop System
CL_gain=dcgain(CL_sys);
% DC-Gain of closed-loop system
K_r=inv(CL_gain); % control signal: u=-K*x+K_r*yd
yd=1; % desired value for output (teta)
%~~~ SIMULATION ~~~
x(:,1)=[0;1;0;0]; % initial condition
for k=2:ct
u(k)=-K*x(:,k-1)+K_r*yd;
x(:,k)=x(:,k-1)+.01*(A*x(:,k-1)+B*u(k));
end
%~~~~~ PLOTTING STATES ~~~~~
figure;subplot(2,2,1);plot(t,x(1,:));
title('Tracking Problem, state X_1');
subplot(2,2,2);plot(t,x(2,:));
title('Tracking Problem, state:X_2');
subplot(2,2,3);plot(t,x(3,:));
title('Tracking Problem, state:X_3');
subplot(2,2,4);plot(t,x(4,:));
title('Tracking Problem, both:X_4 and Output');
%~~~~~ PLOTTING CONTROL SIGNAL ~~~~
figure;plot(t,u);
title('Tracking Problem, Control Signal');