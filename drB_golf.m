function drB = drB_golf(in1,in2)
%drB_golf
%    drB = drB_golf(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    04-Dec-2023 14:38:28

dth1 = in1(3,:);
dth2 = in1(4,:);
l1 = in2(4,:);
l2 = in2(8,:);
th1 = in1(1,:);
th2 = in1(2,:);
t2 = th1+th2;
t3 = cos(t2);
t4 = sin(t2);
drB = [dth1.*(l2.*t3+l1.*cos(th1))+dth2.*l2.*t3;dth1.*(l2.*t4+l1.*sin(th1))+dth2.*l2.*t4;0.0];
end
