function dz = sys(t,z)
%define paraMatrix
A = [0.5444 -0.0020;2.3244 -0.1568];
B = [0;0.8368];
C = [1 2];
D = 0;


el%define setpoint target
if (t>=0 & t<=20)
    z2d = -10;se
    z2d = -20;
end

%define ss function
zd = [0;z2d];
Ke = place(A,B,[-1+j,-1-j]); % Feedback Gain;
e = zd - z;
u = Ke*e;

hold on
plot(t,u);
grid on

dz = A*z + B*u;
end