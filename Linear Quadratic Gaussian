
sys1 = ss(A_lin,[B_lin ones(6,1)],Cx,0); %converting into state space
Qn = (10); % setting covraince matrices
Rn = (0.000001);
[kest,L_kal] = kalman(sys1,Qn,Rn,0);  %using kalman filter to estimate kalman gain
reg = lqgreg(kest, K);
syms x1t(t) x2t(t) x3t(t) x4t(t) x5t(t) x6t(t)
syms x1hat(t) x2hat(t) x3hat(t) x4hat(t) x5hat(t) x6hat(t);
disturbance_vec = [0.001 0.001 0.001 0.001 0.001 0.001]'; % reference and 
control_input3= vpa(-K*(vec_xhat)); % setting control input from LQR
%solving for non linear model by appending states with hat and using ode45
%to solve.
f1_sub4 = simplify(vpa(subs(f1, [F l1 l2 M m1 m2 g],[ control_input3 20 10 1000 100 100 9.8])));
f2_sub4 = simplify(vpa(subs(f2, [F l1 l2 M m1 m2 g],[ control_input3 20 10 1000 100 100 9.8])));
f3_sub4 = simplify(vpa(subs(f3, [F l1 l2 M m1 m2 g],[ control_input3 20 10 1000 100 100 9.8])));
f1_sub4_t = subs(f1_sub4,[ x1 x2 x3 x4 x5 x6],[x1t x2t x3t x4t x5t x6t]);
f2_sub4_t = subs(f2_sub4,[ x1 x2 x3 x4 x5 x6],[x1t x2t x3t x4t x5t x6t]);
f3_sub4_t = subs(f3_sub4,[ x1 x2 x3 x4 x5 x6],[x1t x2t x3t x4t x5t x6t]);
f1_sub4_hat = subs(f1_sub4,[ x1 x2 x3 x4 x5 x6],[x1hat x2hat x3hat x4hat x5hat x6hat]);
f2_sub4_hat = subs(f2_sub4,[ x1 x2 x3 x4 x5 x6],[x1hat x2hat x3hat x4hat x5hat x6hat]);
f3_sub4_hat = subs(f3_sub4,[ x1 x2 x3 x4 x5 x6],[x1hat x2hat x3hat x4hat x5hat x6hat]);
LnC_kal = L_kal*Cx;
diff_eqn11 = diff(x1hat,t,1)==f1_sub4_hat+ LnC_kal(1,:)*(vec_xt-vec_xhat)+disturbance_vec(1,1);
diff_eqn21 = diff(x2hat,t,1)==f2_sub4_hat+ LnC_kal(2,:)*(vec_xt-vec_xhat)+disturbance_vec(2,1);
diff_eqn31 = diff(x3hat,t,1)==f3_sub4_hat+ LnC_kal(3,:)*(vec_xt-vec_xhat)+disturbance_vec(3,1);
diff_eqn41= diff(x4hat,t,1)==x1hat + LnC_kal(4,:)*(vec_xt-vec_xhat)+disturbance_vec(4,1);
diff_eqn51 = diff(x5hat,t,1)==x2hat + LnC_kal(5,:)*(vec_xt-vec_xhat)+disturbance_vec(5,1);
diff_eqn61 = diff(x6hat,t,1)==x3hat+ LnC_kal(6,:)*(vec_xt-vec_xhat)+disturbance_vec(6,1);
diff_eqn71 = diff(x1t,t,1)==  f1_sub4_t + disturbance_vec(1,1);
diff_eqn81 = diff(x2t,t,1)== f2_sub4_t+ disturbance_vec(2,1);
diff_eqn91 = diff(x3t,t,1) == f3_sub4_t+disturbance_vec(3,1);
diff_eqn101= diff(x4t,t,1) == x1t+disturbance_vec(4,1);
diff_eqn111 = diff(x5t,t,1) == x2t+disturbance_vec(5,1);
diff_eqn121 = diff(x6t,t,1) == x3t+disturbance_vec(6,1);
state_vec1= [vec_xhat; vec_xt];
eqns3 = [ diff_eqn11; diff_eqn21; diff_eqn31; diff_eqn41; diff_eqn51; diff_eqn61; diff_eqn71; diff_eqn81; diff_eqn91; diff_eqn101; diff_eqn111;diff_eqn121]; 
[M4,F4] = massMatrixForm(eqns3,state_vec1);
f = M4\F4;
ode_fun2 = odeFunction(f,state_vec1);
x_init2 =[zeros(6,1); x_init];
[t,x_out3] = ode45(ode_fun2,[0 sim_time],x_init2);
figure;
hold on;
%plots
plot(t,x_out3(:,10),'r');
plot(t,x_out3(:,4),'--','Color',[1,0,0]);
plot(t,x_out3(:,11),'g');
plot(t,x_out3(:,5),'--','Color',[0,1,0]);
plot(t,x_out3(:,12),'b');
plot(t,x_out3(:,6),'--','Color',[0,0,1]);
legend('xc','xchat','t1','t1hat','t2','t2hat')
%eig(A_lin-L_kal*Cx) % eigen value check of error
figure;
hold on;
plot(t,x_out3(:,7),'r');
plot(t,x_out3(:,1),'--','Color',[1,0,0]);
plot(t,x_out3(:,8),'g');
plot(t,x_out3(:,2),'--','Color',[0,1,0]);
plot(t,x_out3(:,9),'b');
plot(t,x_out3(:,3),'--','Color',[0,0,1]);
legend('xcdot','xchatdot','t1dot','t1hatdot','t2dot','t2hatdot')

