function sim_golfderive()

    %% Definte fixed paramters
    m1= 0.036; %kg
    m2= 0.0095; %kg
    th1 = -pi/2; %rad
    th2 = 0; %rad
    th1_0 = 0; %rad
    th2_0 = 0; %rad
    dth1 = 0; %rad/s
    dth2 = 0; %rad/s
    l1 = 0.0635; %m
    l2 = 0.0635; %m
    c1 = 0.025; %m
    c2 = 0.034; %m
    g = -10; %m2/s
    I1 = (1/3)*m1*l1^2; %kgm2
    I2 = (1/3)*m2*l2^2; %kgm2
    k = .0994; %Nm

    p   = [m1; I1; c1; l1; m2; I2; c2; l2; g; k; th1_0; th2_0;];       % parameters

    %% Perform Dynamic simulation    
    dt = 0.001;
    tf = 10;
    num_steps = floor(tf/dt);
    tspan = linspace(0, tf, num_steps); 
    z0 = [th1; th2; dth1; dth2];
    z_out = zeros(4,num_steps);
    z_out(:,1) = z0;
    for i=1:num_steps-1
        dz = dynamics(z_out(:,i), p);
        z_out(:,i+1) = z_out(:,i) + dz*dt;
%         z_out(3:4,i+1) = z_out(3:4,i) + dz(3:4)*dt;
%         z_out(1:2,i+1) = z_out(1:2,i) + z_out(3:4,i+1)*dt; % + 0.5*dz(3:4)*dt*dt;
    end
    final_state = z_out(:,end);
    

    %% Compute Energy
    E = energy_golf(z_out,p);
    figure(1); clf
    plot(tspan,E);xlabel('Time (s)'); ylabel('Energy (J)');

    %% Plot theta1 & theta2
    th1 = th1_golf(z_out, p);
    figure(2); clf
    plot(tspan, th1);xlabel('Time (s)'); ylabel('Theta 1');

    th2 = th2_golf(z_out, p);
    figure(3); clf
    plot(tspan, th2);xlabel('Time (s)'); ylabel('Theta 2');

    %% Animate Solution
    figure(4); clf;
        % Prepare plot handles
    hold on
    h_l1 = plot([0],[0],'LineWidth',5);
    h_l2 = plot([0],[0],'LineWidth',3);
    xlabel('x')
    ylabel('y');
    h_title = title('t=0.0s');
    
    axis equal
    axis([-0.3 0.3 -0.3 0.3]);
    skip_frame = 10;
    
    %Step through and update animation
    for i=1:num_steps
        if mod(i, skip_frame)
            continue
        end
        % interpolate to get state at current time.
        t = tspan(i);
        z = z_out(:,i);
        keypoints = keypoints_golf(z,p);

        rB = keypoints(:,1); % Vector to joint 
        rC = keypoints(:,2); % Vector to end effector

        set(h_title,'String',  sprintf('t=%.2f',t) ); % update title
        
        % Plot link 1
        set(h_l1,'XData', [0 ; rB(1)]);
        set(h_l1,'YData', [0 ; rB(2)]);

        % Plot link 2
        set(h_l2,'XData' , [rB(1) ; rC(1)] );
        set(h_l2,'YData' , [rB(2) ; rC(2)] );

        pause(.01)
    end
end

function dz = dynamics(z,p)
    % Get mass matrix
    A = A_golf(z,p);
    
    % Get forces
    u = [0 0]';
    %u = [ 0 -5*z(4)]';
    b = b_golf(z,u, p);
    
    % Solve for qdd
    qdd = A\b;
    dz = 0*z;

    % Form dz
    dz(1:2) = z(3:4);
    dz(3:4) = qdd;
end