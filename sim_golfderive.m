function sim_golfderive()
    k_list = .01:.01:2;
    max_vel = 0;
    max_k = 0;
    for k = k_list
        drB = sim_golfderive_k(.001, false);

        vel_on_impact = sqrt(drB(1)^2 + drB(2)^2 + drB(3)^2);

        if (vel_on_impact > max_vel) 
            max_k = k;
            max_vel = vel_on_impact;
        end
    end

    sim_golfderive_k(max_k, true)

    max_vel
    max_k
    
end

function drB = sim_golfderive_k(k, debug)

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
    %k = .07; %Nm

    p   = [m1; I1; c1; l1; m2; I2; c2; l2; g; k; th1_0; th2_0;];       % parameters

    %% Perform Dynamic simulation    
    dt = 0.00001;
    tf = 0.2;
    num_steps = floor(tf/dt);
    tspan = linspace(0, tf, num_steps); 
    z0 = [th1; th2; dth1; dth2];
    z_out = zeros(4,num_steps);
    z_out(:,1) = z0;
    for i=1:num_steps-1
        dz = dynamics(tspan(i), z_out(:,i), p);
        z_out(:,i+1) = z_out(:,i) + dz*dt;
%         z_out(3:4,i+1) = z_out(3:4,i) + dz(3:4)*dt;
%         z_out(1:2,i+1) = z_out(1:2,i) + z_out(3:4,i+1)*dt; % + 0.5*dz(3:4)*dt*dt;
        theta1 = z_out(1,i);
        theta2 = z_out(2,i);
        thresh = .001;
        %% WANT TO CHANGE THIS SO IT MORE LIKE WHEN IT WOULD HIT THE BALL. ALL GOOD IF IT DOESN"T LINE UP PERFECTLY THAT WILL BE SHOWN IN THE HEAT MAP AS A LOSS.
        %% How about we do position instead of theta values
        current_z = z_out(:,i+1);
        keypoints = keypoints_golf(current_z,p);
        rC = keypoints(:,2); % Vector to end effector
        rC_end_x = rC(1);
        rC_end_y = rC(2);

        % if theta1>0 && (theta2<thresh && theta2 >-thresh) %stops simulation once th1 =0 (pointing down)
        %     t_stop = tspan(i);
        %     num_stop = floor(t_stop/dt);
        %     vx = dth2*(l1+l2);
        %     disp(vx)
        %     break
        % end
        
        % If x value is ~0 and y value is ~-(l1+l2)
        % Seems to work great with just x position
        if (abs(rC_end_x) < thresh)%%(abs(-(l1+l2) - (rC_end_y)) < thresh)
            t_stop = tspan(i);
            num_stop = floor(t_stop/dt);
            vx = dth2*(l1+l2);
            disp(vx);
            break
        end
    end
    final_state = z_out(:,end);

    n = num_stop;
    
    if (debug)
        %% Compute Energy
        E = energy_golf(z_out,p);
        figure(1); clf
        plot(tspan(1:n),E(1:n));xlabel('Time (s)'); ylabel('Energy (J)');
    
        %% Plot theta1 & theta2
        th1 = th1_golf(z_out, p);
        figure(2); clf
        plot(tspan(1:n), th1(1:n));xlabel('Time (s)'); ylabel('Theta 1');
    
        th2 = th2_golf(z_out, p);
        figure(3); clf
        plot(tspan(1:n), th2(1:n));xlabel('Time (s)'); ylabel('Theta 2');
    
        %% Plot velocities
    
        dth1 = z_out(3, 1:num_stop);
        dth2 = z_out(4, 1:num_stop);
    
        figure(5); clf
        plot(tspan(1:num_stop), dth1, 'b-');
        xlabel("Time (s)"); ylabel("dth1 (Rad/S)");
        figure(6); clf
        plot(tspan(1:num_stop), dth2, 'k-');
        xlabel("Time (s)"); ylabel("dth2 (Rad/S)");
    
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
        for i=1:num_stop
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
    
    z_out = z_out(:, 1:num_stop);

    %% How do we get back just the velocity of the end-point
    drB = drB_golf(z_out(:,num_stop), p);
end

function tau = control_law(t, z, p)
    tau1_des = torque_curve(z(3)); %[Nm] desired torque to be applied
    tau1_t = 0.01; %[s] time torque will be applied

    %stall torque - torque constant * speed

    if t < tau1_t
        tau = [tau1_des 0]';
    else
        tau = [0 0]';
    end
    
end

function dz = dynamics(t,z,p)
    % Get mass matrix
    A = A_golf(z,p);

    % Compute Controls 
    tau = control_law(t,z,p);
    
    % Get b = Q - V(q,qd) - G(q)
    b = b_golf(z,tau,p);
    
    % Solve for qdd
    qdd = A\b;
    dz = 0*z;

    % Form dz
    dz(1:2) = z(3:4);
    dz(3:4) = qdd;
end

function tau_m = torque_curve(dth1_input)
    % dth1_input should be rad/s.
    % first convert to rpm
    rpm = dth1_input * (60)/(2*pi);

    tau_m = 320/2.2 - (1/2.2)*rpm;
end