%close all
time_step = 0.01
sim_time = 100
syms M m1 m2 l1 l2 x1 x2 x3 x4 x5 x6 g
syms F real
syms M m1 m2 l1 l2 x1 x2 x3 x4 x5 x6 g real
syms k1 k2 k3 k4 k5 k6 real

%Representing the non-linear state-space equations
f1 = (F - 0.5* m1*g*sin(2*x5) - m1*l1*sin(x5)*x2^2-0.5*m2*g*sin(2*x6)-m2*l2*sin(x6)*x3^2)/(M+m1*sin(x5)^2+m2*sin(x6)^2);
f2 = ((F - m1*l1*sin(x5)*x2^2 - m2*l2*sin(x6)*x3^2 - m2*g*sin(x6)*cos(x6))*cos(x5) - (M+m1+m2*sin(x6)^2)*g*sin(x5))/(l1*(M+m1*sin(x5)^2+m2*sin(x6)^2));
f3 = ((F - m1*l1*sin(x5)*x2^2 - m2*l2*sin(x6)*x3^2 - m1*g*sin(x5)*cos(x5))*cos(x6) - (M+m2+m1*sin(x5)^2)*g*sin(x6))/(l2*(M+m1*sin(x5)^2+m2*sin(x6)^2));
f4 = x1;
f5 = x2;
f6 = x3;
K = [ k1 ;k2; k3; k4; k5; k6];

%Linearizing the state-dynamics matrix A, B

%Linearize the matrix A at the equilibrium point
Jacob_A = [ diff(f1,x1) diff(f1,x2) diff(f1,x3) diff(f1,x4) diff(f1,x5) diff(f1,x6);
    diff(f2,x1) diff(f2,x2) diff(f2,x3) diff(f2,x4) diff(f2,x5) diff(f2,x6);
    diff(f3,x1) diff(f3,x2) diff(f3,x3) diff(f3,x4) diff(f3,x5) diff(f3,x6);
    diff(f4,x1) diff(f4,x2) diff(f4,x3) diff(f4,x4) diff(f4,x5) diff(f4,x6);
    diff(f5,x1) diff(f5,x2) diff(f5,x3) diff(f5,x4) diff(f5,x5) diff(f5,x6);
    diff(f6,x1) diff(f6,x2) diff(f6,x3) diff(f6,x4) diff(f6,x5) diff(f6,x6) ];
%Substituting the values of equilibrium point
Jacob_A_sub=subs(Jacob_A, [x1 x2 x3 x4 x5 x6],[ 0 0 0 0 0 0]);
D = M + m1*sin(x5)^2 + m2*sin(x6)^2;


%Linearize the matrix B at the equilibrium point
Jacob_B = [ diff(f1,F); diff(f2,F); diff(f3,F); diff(f4,F); diff(f5,F); diff(f6,F)];
%Substituting the values of equilibrium point
Jacob_B_sub = subs(Jacob_B, [x1 x2 x3 x4 x5 x6],[ 0 0 0 0 0 0]);
closed_loop = Jacob_B_sub*K'
controllability_test = [ Jacob_B_sub Jacob_A_sub*Jacob_B_sub Jacob_A_sub^2*Jacob_B_sub Jacob_A_sub^3*Jacob_B_sub Jacob_A_sub^4*Jacob_B_sub Jacob_A_sub^5*Jacob_B_sub];
controllability_test_sub = double(subs(controllability_test,[l1 l2 M m1 m2 g],[20 10 1000 100 100 10]));
check_rank = rank(controllability_test_sub)

%
A_lin1 = subs(Jacob_A, [l1 l2 M m1 m2 g],[20 10 1000 100 100 9.8]);
A_lin = double(subs(A_lin1,[x1 x2 x3 x4 x5 x6],[0 0 0 0 0 0]));
B_lin1 = subs(Jacob_B,[l1 l2 M m1 m2 g ],[20 10 1000 100 100 9.8]);
B_lin = double(subs(B_lin1,[x1 x2 x3 x4 x5 x6],[0 0 0 0 0 0]));

% Q_cost = 0.1*[1 5 5 5 5 5
%           5 1 5 5 5 5 
%           5 5 1 5 5 5
%           5 5 5 1 5 5
%           5 5 5 5 1 5
%           5 5 5 5 5 1];
%Q_cost = 1000*eye(6,6);

%  Q_cost = 10*[10 3 2 5 2 3;
%               3 100 5 2 2 3;
%               2 4 5 2 2 4;
%               5 7 6 5 6 6;
%               2 5 2 5 100 6;
%               3 7 3 5 4 10]; 
%   R_cost = 0.0001;          %Q1_R1

% Q_cost = 10*[20 0 0 0 0 0;
%               0 100 0 0 0 0;
%               0 0 10 0 0 0;
%               0 0 0 10 0 0;
%               0 0 0 0 100 0;
%               0 0 0 0 0 20];
% 
% R_cost = 0.00002;         %Q2_R2

% Q_cost = 10*[1 0 0 0 0 0;
%                0 1 0 0 0 0;
%                0 0 1 0 0 0;
%                0 0 0 1 0 0;
%                0 0 0 0 1 0;
%                0 0 0 0 0 1];
% R_cost = 1;               %Q3_R3
% Q_cost = 1*[0.5 0 0 0 0 0 ;
%              0 0.5 0 0 0 0;
%              0 0 0.5 0 0 0;
%              0 0 0 0.5 0 0;
%              0 0 0 0 0.5 0;
%              0 0 0 0 0 0.5];
% R_cost = 0.5;             %Q4_R4

%Defining the positive-definite weight matrices Q and R
Q_cost = 10*[10 3 2 5 2 3;
              3 100 5 2 2 3;
              2 4 5 2 2 4;
              5 7 6 5 6 6;
              2 5 2 5 100 6;
              3 7 3 5 4 10]; 
R_cost = 0.0001;          %Q5_R5

% Q_cost = 1*[10 5 80 2 6 10 ;
%              10 20 100 50 20 10;
%              30 50 30 10 10 10;
%              20 60 10 30 20 40;
%             30 10 10 10 20 80;
%              60 10 20 60 50 10];
% R_cost = 0.001;     


%Designing the LQR controller, by defining the Q and R matrices above,
%and giving A,B,Q,R as the parameters of lqr command.

[K,s,e]= lqr(A_lin,B_lin,Q_cost,R_cost,0);
syms t real
x_init = 0.001*ones(6,1);
x_fin=[];
test= A_lin-B_lin*K;
eig(test)
for t=0:time_step:sim_time
    x_fin_temp = expm((A_lin-B_lin*K).*t)*x_init;
    x_fin= [x_fin x_fin_temp];
end

%figures for states-position, theta1 and theta2
figure
title('States of the system')
hold on
plot([0:time_step:sim_time],x_fin(1,:),'b');
plot([0:time_step:sim_time],x_fin(2,:),'g');
plot([0:time_step:sim_time],x_fin(3,:),'r');
legend('Cart-position(linear)','Pendulum1(linear)','Pendulum2(linear)')



% Defining the control input 
u=K*x_fin;
figure
hold on
title('Control Input')
plot(u,'b')






