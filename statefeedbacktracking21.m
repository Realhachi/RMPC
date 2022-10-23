A = [0.5444 -0.0020;2.3244 -0.1568];
B = [0;0.8368];
C = [0 1];
D = 0;
P=[-15+10i -15-10i];
% Desired Closed-Loop Ploes
K=acker(A,B,P); % State-Feedback Gain
sys = ss(A-B*K,B,C,D);
t=0:0.01:10;
x = initial(sys,[1;10],t);
x1 = [1 0]*x';
