function keypoints = keypoints_golf(in1,in2)
%KEYPOINTS_GOLF
%    KEYPOINTS = KEYPOINTS_GOLF(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    15-Nov-2023 14:55:17

c1 = in2(3,:);
c2 = in2(7,:);
l1 = in2(4,:);
l2 = in2(8,:);
th1 = in1(1,:);
th2 = in1(2,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = l1.*t2;
t6 = l1.*t3;
t7 = -t5;
keypoints = reshape([t6,t7,t6+l2.*sin(th2),t7-l2.*cos(th2),c1.*t3,-c1.*t2,t6+c2.*sin(t4),t7-c2.*cos(t4)],[2,4]);
end