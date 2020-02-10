close all
clear variables
clc


gravity = [-8.00665 0 0]';
mass_dry = 27215; 
mass_wet = 38963; 

Isp = 282;

thrbk_r = 0.8;
% rho_1 = 400000;
rho_2 = 845200; 
rho_1 = rho_2 * thrbk_r; 
position_0 = [5582 2100 450]'; 
velocity_0 = [-353 -190 -40]';

gamma = pi*(40/180);

% %optimal final time
% tf_opt = 30;
% %discrete time steps
% N = tf_opt*2;
% %time increment
% time_step = tf_opt/N;

%discrete time steps
N = 55;
%optimal final time
tf_opt = 30;
%time increment
time_step = tf_opt/N;

%allowed landing error(3m and 0.5m/s of horizontal speed)
land_error = 3;
land_error_v = 0.5;

alph = -1/gravity(1)/Isp;

%Convex Linear Operators
E = [eye(3) zeros(3,4)];
F = [0 0 0 0 0 0 1];
E_u = [eye(3) zeros(3,1)];
E_v = [zeros(3,3) eye(3) zeros(3,1)];

e_1 = [1 0 0]';
e_2 = [0 1 0]';
e_3 = [0 0 1]';

S = [e_2';e_3'];

c = e_1/tan(gamma);

%Define the State transition matrices 
A_c = [zeros(3,3) eye(3) zeros(3,1);zeros(4,7)];
B_c = [[zeros(3,3);eye(3);0 0 0] [0 0 0 0 0 0 -alph]'];

A = expm(A_c*time_step); %continuous time A matrix

B = A*(eye(7)*time_step - A_c*time_step^2/2)*B_c; %continuous time B matrix

Lambda_k = zeros(N*7,4);
Lambda_k(1:7,1:4) = B;
for k = 2:N
    
    Lambda_k((k-1)*7+1:k*7,:) = A*Lambda_k((k-2)*7+1:(k-1)*7,:) + B; 
end

Psi_k = zeros(N*7,4*(N+1));
for k = 2:(N)
    %Next time step of gravities effect is dependent on the state
    %transition matrix and the previous time step
    Psi_k((k-1)*7+1:k*7,:) = A*Psi_k((k-2)*7+1:(k-1)*7,:);
    Psi_k((k-1)*7+1:k*7,((k*4-7):(4*k-4))) = B;
end

% Mass after the change of variables
z0 = log(mass_wet-alph*rho_2*time_step*(0:N)');

% Initial state vector
y0 = [position_0; velocity_0; log(mass_wet)];

s(1:N,7) = 0;
for i = 1:N
    s(i,:) = (7*i-6):(7*i);
end

cvx_begin
    variable eta((N+1)*4)
    variable y(N*7)

    % Objective function minimun error
%     minimize(norm(y(end-6:end-4),2))
    % Objective function minimun fuel consumption
    minimize(y(N*7))

    subject to

    % Convexified thrust constraint
    for k = 0:N
        norm(E_u*eta(4*k+1:4*k+4), 2) <= eta(4*k+4);
    end

    % Thrust constraint 1
    eta(4) <= rho_2*exp(-z0(1)).*(1-(F*y0-z0(1)));
    rho_1*exp(-z0(1))*(1-(F*y0-z0(1))+0.5*(F*y0-z0(1)).^2) <= eta(4);

    for k = 1:N
        % Cone constraints
        norm(S*E*(y(s(k, :))-y(s(N, :))), 2)-c'*(E*(y(s(k, :)))) <= 0;

        % Thrust constraints
        eta(4*(k)+4) <= rho_2*exp(-z0(k+1)).*(1-(F*y(s(k, :))-z0(k+1)));

        rho_1*exp(-z0(k+1))*(1-(F*y(s(k, :))-z0(k+1))+...
            0.5*(F*y(s(k, :))-z0(k+1)).^2) <= eta(4*k+4);

        % System dynamics constraint
        y(s(k, :)) == A^k*y0+Lambda_k(s(k, :), :)*[gravity; 0]+...
                        Psi_k(s(k, :), :)*eta;
     end

    % Fuel mass constraint
    y(end) >= log(mass_dry);

    % Final height is 0 constraint
    y(end-6) == 0;

    % Final velocity constraint
    for i = 1:3
        y(end-i) == 0;
    end
    
    % Trajectory constraint
    for i = 0:(N-1)
        % No -y -z
        y(7*i+2) >= 0;
        y(7*i+3) >= 0;
        % No go up
        y(7*i+4) <= 0;
    end
    
    % Final velocity & position
    norm(y(end-6:end-4),2) <= land_error;
    norm(y(end-9:end-8),2) <= land_error_v;
cvx_end

% Converting output into manageable format
dist(1:3, N+1) = 0;
dist(1:3, 1) = position_0;
mass(1) = mass_wet;
for i = 1:N
    dist(1:3, i+1) = y((7*i-6):(7*i-4));
    mass(i+1) = y(7*i);
end

acc(N+1, 4) = 0;
velocity(N+1, 4) = 0;
for i = 1:N
    acc(i, 1:3) = -eta(i*4-3:i*4-1);
    acc(i, 4) = -norm(acc(i, 1:3), 2);
    velocity(i, 1:3) = y(i*7-3:i*7-1);
    velocity(i, 4) = norm(velocity(i, 1:3), 2);
end

thrust(N+1, 4) = 0;
thrust(1, :) = -mass(1).*acc(1, :);
for i = 1:N
    thrust(1+i, :) = -exp(mass(1+i)).*acc(1+i, :);
end

figure
subplot(2,3,1)
plot(0:time_step:tf_opt,dist)
grid on
legend('x','y','z');
xlabel('Time (seconds)')
ylabel('Position (meters)')
title('Trajectory vs. Time')
subplot(2,3,2)
plot(0:time_step:tf_opt,velocity(:,1:3))
grid on
xlabel('Time (seconds)')
ylabel('velocity (m/s)')
title('velocity vs. Time')
subplot(2,3,4)
plot(0:time_step:tf_opt, thrust(:, 1:3))
grid on
xlabel('Time (seconds)')
ylabel('Thrust/N')
title('thrust vs. Time')
grid on
subplot(2,3,5)
plot(0:time_step:tf_opt, thrust(:, 4))
grid on
xlabel('Time (seconds)')
ylabel('Thrust/N')
title('Thrust vs. Time')
grid on

subplot(2,3,[3 6])
plot3(dist(3, :), dist(2, :), dist(1, :), 'r.')
grid on
xlabel('z (meters)')
ylabel('y (meters)')
zlabel('x (meters)')
title('3D Trajectory')
axis equal
