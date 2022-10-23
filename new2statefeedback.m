clear 
close all

m = 4; J = 0.0475; r = 0.25; g = 9.8; c = 0.05; F1=0; F2 = m*g;
A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 0 -g -c/m 0 0;
     0 0 0 0 -c/m 0;
     0 0 0 0 0 0];
B = [0 0;
      0 0;
      0 0;
      1/m 0;
      0 1/m;
      r/J 0];
C = [1 0 0 0 0 0;
       0 1 0 0 0 0];
Qz = eye(6); Qv = eye(2);
 
figure()
sys = ss(A, B, C, 0);
step(sys, 10)


%%%%%%%%%%%%%%%% Verbose mode: how to plot output response %%%%%%%%%%%%%%%%
syms s
x0 = [0.01; 0.01; 0; 0; 0; 0];   % initial condition
X_zs = (s*eye(6) - A)\B* (1/s);  % Unit step response.  This is called zerso state response. The matrix (sI - A)^-1 is the state transition matrix, e^(A*t)
X_zi = (s*eye(6) - A)\x0;       % This is zero input response
out = C * ( ilaplace(X_zs) + ilaplace(X_zi) ); % since system is linearized, output = X_zs + X_zi (super position principle)

figure()
fplot(out(1,1), [0 10])  % plot for 10 sec
figure()
fplot(out(2,2), [0 10])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%  Trajectory tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = 0:0.1:10;
 
figure()
xd = [1; 0; 0; 0; 0; 0]; yd = [0; 1; 0; 0; 0; 0]; % desired values of x, y
K = lqr(A, B, Qz, Qv);    % linear quadratic regulator
H1ax = ss(A-B*K,B(:,1)*K(1,:)*xd,C(1, :),0);  % u = -K (x - xd); x_dot = (A - B*K)x + B*K*xd; y = C*x
H1ay = ss(A-B*K,B(:,2)*K(2,:)*yd,C(2, :),0);
step(H1ax, H1ay, 10);
legend('x', 'y');