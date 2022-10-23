% LQ Tracking Control
clear; close all; clc;

velocity_result=load('velocity_result.txt');

h=0.01;
M=6762;

A=1;
B=h/M;
C=1;
% A=[0,1;0,0];
% B=[0;1/M];
% C=[0,1];

S=1;
%Q=10000000;
Q=100000000;
R=0.0001;

% WP=load('Path_rotated.txt');
%Tf=400;
[final_index,~]=size(velocity_result(:,1));
v=zeros(final_index,1);
v(1)=0;

%t=0:T:50;
%t=velocity_result(:,1);
time=zeros(final_index,1);
velocity_by_time=velocity_result(:,3);
% y(1:100)=10;
% y(101:200)=10;
% y(201:501)=1;
% y=10*sin(t);
u=zeros(final_index-1,1);
iteration=40;        % number of iteration

for j=1:final_index-iteration-1
    fprintf('j = %d / %d\n',j, final_index);
    input_velocity=velocity_by_time(j:j+iteration-1);
    
    H=C'*S*C;
    g=-C'*S*input_velocity(iteration);
    
    for k=iteration-1:-1:1       % Riccati equation
        g=(A-B*inv(R+B'*H*B)*B'*H*A)'*g-C'*Q*input_velocity(k);
        H=A'*H*A-A'*H*B*inv(R+B'*H*B)*B'*H*A+C'*Q*C;
    end
    u(j)=-inv(R+B'*H*B)*B'*(H*A*v(j)+g);
    
    v(j+1)=A*v(j)+B*u(j);   % ´ÙÀ½ step¿¡¼­ÀÇ ¼Óµµ
    time(j+1)=j*h;
    
end
subplot(2,1,1)
plot(time(1:final_index-iteration),velocity_by_time(1:final_index-iteration),'--b')
hold on
plot(time(1:final_index-iteration),v(1:final_index-iteration),'r')
xlabel('time[s]')
ylabel('m/s')
legend('v_d','v_c')
grid on
hold off

subplot(2,1,2)
plot(time(1:final_index-iteration-1),u(1:final_index-iteration-1))
xlabel('time[s]')
ylabel('F_x_d [N]')
grid on
hold off