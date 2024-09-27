clear;
clc;
% MOI
global J R_rw_to_body;
J= [100 30 -10;
    30 150 20;
    -10 20 170]; %Total inertia
dt=0.1;
T=60;
time=0:dt:T;
n=length(time);
qd= [0.640381703895931 -0.312136686798122 0.701770590770257 0 ]; % desired attitude
qcm=[qd(4) qd(3) -qd(2) qd(1);
    -qd(3) qd(4) qd(1) qd(2);
     qd(2) -qd(1) qd(4) qd(3);
    -qd(1) -qd(2) -qd(3) qd(4)];
q0=[0 0 0 1];
T_rw(:,1)=[0 0 0]'; %torque at time zero
qe= (qcm*[-q0(1);-q0(2);-q0(3);q0(end)])'; % qe=qcm*q0^-1
x=zeros(10,n); % quaternion; omega;H_rw 
x(4,1)=1; % initial attitude the same problem 1~5.
R_rw_to_body=eye(3,3); % RW configuration along the principal axis for problem 1~4.
%% Problem 1
kp= (6.6/9.1)^2; kd=13.2/9.1;
for i = 2:n
    T_comm= 2.*kp.*J*qe(i-1,1:3)'-kd.*J*x(5:7,i-1);
    T_rw(:,i)= R_rw_to_body' * T_comm;
    x(:,i)=RK4(@xdot, x(:,i-1),dt,T_rw(:,i));
    qe(i,:)= (qcm*[-x(1,i);-x(2,i);-x(3,i);x(4,i)])';
end
mustplot_P2(time, x(1:4,:), x(5:7,:), T_rw, x(8:end,:), R_rw_to_body');
%% Problem 2
x(5:7,1)=[deg2rad(1) deg2rad(1) deg2rad(1)]'; %initial Angular Vel.
for i = 2:n
    T_comm= 2.*kp.*J*qe(i-1,1:3)'-kd.*J*x(5:7,i-1);
    T_rw(:,i)= R_rw_to_body' * T_comm;
    x(:,i)=RK4(@xdot, x(:,i-1),dt,T_rw(:,i));
    qe(i,:)= (qcm*[-x(1,i);-x(2,i);-x(3,i);x(4,i)])';
end
mustplot_P2(time, x(1:4,:), x(5:7,:), T_rw, x(8:end,:), R_rw_to_body');
%% Problem 3
x(5:7,1)=[0 0 0]';
kp=(6.705/23)^2; kd=13.2/22.9;
for i = 2:n
    T_comm= 2.*kp.*J*qe(i-1,1:3)'-kd.*J*x(5:7,i-1);
    temp_T_rw=R_rw_to_body' * T_comm; % temporary T_rw
    if (norm(temp_T_rw)<2)
         T_rw(:,i)=temp_T_rw;
    else
        T_rw(:,i)= 2/norm(temp_T_rw) *temp_T_rw;
    end
    x(:,i)=RK4(@xdot, x(:,i-1),dt,T_rw(:,i));
    qe(i,:)= (qcm*[-x(1,i);-x(2,i);-x(3,i);x(4,i)])';
end
mustplot_P2(time, x(1:4,:), x(5:7,:), T_rw, x(8:end,:), R_rw_to_body');
%% Problem 4
kp=(6.6/19.25)^2; kd=13.2/19;
for i = 2:n
    T_comm= 2.*kp.*J*qe(i-1,1:3)'-kd.*J*x(5:7,i-1);
    temp_T_rw=R_rw_to_body' * T_comm;
    m= 2/max(abs(temp_T_rw));
    T_rw(:,i)=m*temp_T_rw;
    x(:,i)=RK4(@xdot, x(:,i-1),dt,T_rw(:,i));
    qe(i,:)= (qcm*[-x(1,i);-x(2,i);-x(3,i);x(4,i)])';
end
mustplot_P2(time, x(1:4,:), x(5:7,:), T_rw, x(8:end,:), R_rw_to_body');
%% Problem 5
% *** First Method: Maximizing the momentum along axis of rotation ****
% z=[0.6404 -0.3129 0.7018];
% y=[0.5 0.5 -0.233]; y=y/norm(y);
% X=cross(y,z); 
% T_body_to_intermediate=[X; y; z];
% T_intermediate_to_rw=quaternion_to_dcm([sind(56.6/2)/sqrt(2) -sind(56.6/2)/sqrt(2) 0 cosd(56.6/2)]);
%R_rw_to_body=(T_intermediate_to_rw' )*(T_body_to_intermediate');
% *** ends here***
%  ****** Second Method: Randomly putting RWA and Choosing the best******
%R_rw_to_body=euler_to_dcm([rand(1)*360-180,rand(1)*180-90,rand(1)*360-90])
R_rw_to_body=[-0.259312660416426   0.064401729200469   0.963643793849026;
   0.811476306280041               0.555571741774174   0.181235327925082;
  -0.523701392487987               0.828970721447989  -0.196327263741292];
kp=(6.6/16.72)^2; kd=13.351/16.45;
%  ****** ends here ********
for i = 2:n
    T_comm= 2.*kp.*J*qe(i-1,1:3)'-kd.*J*x(5:7,i-1);
    temp_T_rw=R_rw_to_body' * T_comm;
    m= 2/max(abs(temp_T_rw));
    T_rw(:,i)=m*temp_T_rw;
    x(:,i)=RK4(@xdot, x(:,i-1),dt,T_rw(:,i));
    qe(i,:)= (qcm*[-x(1,i);-x(2,i);-x(3,i);x(4,i)])';
end
mustplot_P2(time, x(1:4,:), x(5:7,:), T_rw, x(8:end,:), R_rw_to_body');
%% Integrator
function xo = RK4(Func, x, dt, u)
    % u is input, x is state
    k1  = Func(x, u)*dt;
    k2  = Func(x + k1*0.5, u)*dt;
    k3  = Func(x + k2*0.5, u)*dt;
    k4  = Func(x + k3 , u)*dt;

    xo = x + (k1 + 2*k2 + 2*k3 + k4)/6;
end
function out = xdot(x, T_rw) 
    global J R_rw_to_body;
    q = x(1:4,1);
    w = x(5:7,1);
    H_rw = x(8:end,1);
    out(1:4,1)=quater(q,w);
    T_comm=R_rw_to_body*T_rw;
    out(5:7,1)= inv(J)* (T_comm-cross(w,J*w+R_rw_to_body*H_rw));
    out(8:10,1)=-T_rw;
end
function [qdot] = quater(q,w)
Omega = [0 w(3) -w(2) w(1);-w(3) 0 w(1) w(2); w(2) -w(1) 0 w(3); -w(1) -w(2) -w(3) 0];
qdot = 1/2*Omega * q;
end