name = 'golf';

% Define variables for time, generalized coordinates + derivatives, controls, and parameters 
syms t g th1 th2 dth1 dth2 ddth1 ddth2 m1 I1 c1 l1 m2 I2 c2 l2 tau1 tau2 k th1_0 th2_0 real

% Group them
q   = [th1; th2];      % generalized coordinates
dq  = [dth1; dth2];    % first time derivatives
ddq = [ddth1; ddth2];  % second time derivatives
u   = [tau1; tau2];     % controls
p   = [m1; I1; c1; l1; m2; I2; c2; l2; g; k; th1_0; th2_0;];        % parameters

% Generate Vectors and Derivatives
ihat = [1; 0; 0];
jhat = [0; 1; 0];
khat = cross(ihat,jhat);

ddt = @(r) jacobian(r,[q;dq])*[dq;ddq]; % a handy anonymous function for taking time derivatives

rA = l1*sin(th1)*ihat + l1*cos(th1)*(-jhat);
rB = rA + l2*sin(th1+th2)*ihat + l2*cos(th1+th2)*(-jhat);
rc1 = c1*sin(th1)*ihat + c1*cos(th1)*(-jhat);
rc2 = rA + c2*sin(th1+th2)*ihat + c2*cos(th1+th2)*(-jhat);

drA = ddt(rA);
drB = ddt(rB);
drc1 = ddt(rc1);
drc2 = ddt(rc2);

% Calculate Kinetic Energy, Potential Energy, and Generalized Forces
F2Q = @(F,r) simplify(jacobian(r,q)'*(F));    % force contributions to generalized forces
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));   % moment contributions to generalized forces

T1 = (1/2)*m1*dot(drc1, drc1) + (1/2)*I1*dth1^2;
T2 = (1/2)*m2*dot(drc2, drc2) + (1/2)*I2*dth2^2;

Vg1 = m1*g*dot(rc1, -jhat);
Vg2 = m2*g*dot(rc2, -jhat);
Ve1 = (1/2)*k*(th2)^2;

T = simplify(T1 + T2);
V = simplify(Vg1 + Vg2 + Ve1);

Q1 = M2Q(tau1*khat, dth1*khat);
Q2 = M2Q(tau2*khat, dth2*khat);
Q = Q1+Q2;

% Assemble the array of cartesian coordinates of the key points
keypoints = [rA(1:2) rB(1:2) rc1(1:2) rc2(1:2)];

%% All the work is done!  Just turn the crank...
% Derive Energy Function and Equations of Motion
E = T+V;
L = T-V;
g = ddt(jacobian(L,dq).') - jacobian(L,q).' - Q;


% Rearrange Equations of Motion
A = jacobian(g,ddq);
b = A*ddq - g;

% Write Energy Function and Equations of Motion
z  = [q ; dq];
matlabFunction(A,'file',['A_' name],'vars',{z p});
matlabFunction(b,'file',['b_' name],'vars',{z u p});
matlabFunction(E,'file',['energy_' name],'vars',{z p});
matlabFunction(keypoints,'file',['keypoints_' name],'vars',{z p});
matlabFunction(th1,'file',['th1_' name],'vars', {z p});
matlabFunction(th2, 'file', ['th2_' name],'vars', {z p});