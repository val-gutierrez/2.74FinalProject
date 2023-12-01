function sim_golfderive()
    % k_list = 0:.001:.5;
    % tau1_t_list = .1;
    % max_vel = 0;
    % max_k = 0;
    % vels = [];
    % mE = 0;
    % for k = k_list
    %     drB = sim_golfderive_k(k, mE, false);
    % 
    %     vel_on_impact = drB(1);%sqrt(drB(1)^2 + drB(2)^2 + drB(3)^2);
    % 
    %     vels = [vels, vel_on_impact];
    % 
    %     if (vel_on_impact > max_vel) 
    %         max_k = k;
    %         max_vel = vel_on_impact;
    %     end
    % end
    % 
    % sim_golfderive_k(max_k, mE, true)
    % vels
    % 
    % max_vel
    % max_k
    % figure(7); clf
    % plot(k_list, vels);
    % xlabel('k value'); ylabel('Impact Velocity');


    %%% NEW CODE AS OF 11/30/23 ----------------------

    % Resonably can do 0 - 1.5
    k_list = 0:.005:.5;
    % Resonably can do 0 - .1
    mE_list = 0:.01:.1;
    % Resonably can do 0 - .2
    tau1_t_list = .2:-.01:0; 

    max_vel = 0;
    opt_k = 0;
    opt_mE = 0;
    opt_tau1_t = 0;

    heat_map_values = zeros(length(k_list), length(mE_list));

    debug = false;

    for k = k_list
        for mE = mE_list
            % best veloctiy and tau_t value in range for K and mE
            max_tau_vel = 0;
            max_tau = 0;
            % Need to run for several tau_1's and get the best value.
            % Discard ones that don't work
            for tau1_t = tau1_t_list
                drB = sim_golfderive_k(k, mE, tau1_t, debug);
                
                % If drB(1) is -1 or the conditions weren't met with tau1_t
                % to get valid results then ignore that point and move on
                % to the next
                if(drB(1) == -1)
                    % Ignore data
                else
                    % data is good, now check if best
                    if (drB(1) > max_tau_vel)
                        max_tau_vel = drB(1);
                        max_tau = tau1_t;
                    end
                end
            end

            % Done running tau1_t list now get best value and store in
            % matrix
            i = find(k_list == k);
            j = find(mE_list == mE);
            heat_map_values(i, j) = max_tau_vel;

            % Checking for best optimization for all
            if(max_tau_vel > max_vel)
                max_vel = max_tau_vel;
                opt_k = k;
                opt_mE = mE;
                opt_tau1_t = max_tau;
            end
        end
    end

    %% Once all data points are gotten

    % 1 lets plot our heat map
    figure(7); clf
    surf(k_list, mE_list, heat_map_values');
    xlabel('k');
    ylabel('mE');
    zlabel('Vel');
    title('3D Surface Plot of Impact Velocity');

    % 2 print best values
    max_vel
    opt_k
    opt_mE
    opt_tau1_t
end

function drB = sim_golfderive_k(k, mE, tau1_t, debug)

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

    p   = [m1; I1; c1; l1; m2; I2; c2; l2; mE; g; k; th1_0; th2_0;];       % parameters

    %% Perform Dynamic simulation  
    num_stop = -1;
    dt = 0.00001;
    tf = 0.3;
    num_steps = floor(tf/dt);
    tspan = linspace(0, tf, num_steps); 
    z0 = [th1; th2; dth1; dth2];
    z_out = zeros(4,num_steps);
    z_out(:,1) = z0;
    for i=1:num_steps-1
        dz = dynamics(tspan(i), z_out(:,i), p, tau1_t);
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
        
        % If x value is ~0 and theta value is ~0 Stop
        if (abs(rC_end_x) < thresh) && (abs(theta1) < thresh)
            t_stop = tspan(i);
            num_stop = floor(t_stop/dt);
            vx = dth2*(l1+l2);
            disp(vx);
            break
        end
    end
    
    if (num_stop == -1)
        %% Discard data for that run
        drB = -1;
        n = num_steps;
    else
        %% Collect data
        final_state = z_out(:,num_stop);
        
        % if final velocity is greater than last recorded update
        drB = drB_golf(final_state, p);

        n = num_stop;
    end
    
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
end

function tau = control_law(t, z, tau1_t)
    tau1_des = torque_curve(z(3)); %[Nm] desired torque to be applied
    % tau1_des = torque_curve(z(3)); %[Nm] desired torque to be applied
    % tau1_t = 0.01; %[s] time torque will be applied

    %stall torque - torque constant * speed

    if t < tau1_t
        tau = [tau1_des 0]';
    else
        tau = [0 0]';
    end
    
end

function dz = dynamics(t,z,p,tau1_t)
    % Get mass matrix
    A = A_golf(z,p);

    % Compute Controls 
    tau = control_law(t,z,tau1_t);
    
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
    %rpm = dth1_input * (60)/(2*pi);

    tau_m = ((85/1000)*(1/9.8)) - ((85/1000)*(1/9.8))/(504/(2*pi)*60)*dth1_input;
end