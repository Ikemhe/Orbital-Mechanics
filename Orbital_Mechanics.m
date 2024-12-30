% Serena I. Elijah
% AE 313 - Space Mechanic HW # 1 
%Professor Luis Mendoza
clc; clear; close all;

%Given Initial Conditions

r_0 = [8000;0;6000] ; %km
v_0 = [0;8.5;0] ; %km/s
t_initial = 0; % Initial Time [s]
t_step = 10; % Time Step [s]


% Constants
Re = 6371; % Mean Radius of the Earth [km] 

%A)1.

h_0 = cross(r_0,v_0);  %Orbit's angular momentum

G_0 = 6.674 * 10^(-20) ;  %Gravitational constant
 
me = 5.9722*(10^24) ; %mass of the Earth in kg
mui = G_0*me ; %we do not ass the mass of the spacecraft because it is negligible compared to the Earth's mass
sr = r_0.*r_0 ;
dr = sum(sr) ;
magnr_0 = sqrt(dr) ; %magnitude of r_0

%A) 2.
e_0 = ((cross(v_0,h_0))/mui) - (r_0/magnr_0) ;

sh = h_0.*h_0 ;
dh = sum(sh) ;
magnh_0 = sqrt(dh) ; %magnitude of h_0

se = e_0.*e_0 ;
de = sum(se) ;
magne_0 = sqrt(de) ; %magnitude of e_0

%A)3.
r_P = ((magnh_0)^2/mui)/(1+(magne_0)) ;%Distance from Earth's center at perigee

%A)4.
r_A = ((magnh_0)^2/mui)/(1-(magne_0)) ; %Distance from Earth's center at apogee

%A)5.
a = (r_A + r_P)/2 ; %semimajor axis

%A)6.
T = (2*pi)*(mui^(-1/2))*(a^(3/2)) ; %s %Orbital period

t_final = T; % Final Time [s]

%how to display one row vectors
%x = [1, 2, 3]
%disp(sprintf('Answer: (%d, %d, %d)', x))


fprintf('The obatined Orbital angular momentum vector is %0.2f,%0.2f,%0.2f\n',h_0)
fprintf('The obatined Orbital eccentricity vector is %0.4f,%0.2f,%0.4f\nThe magnitude of the eccentricity vector is %0.4f, and it is <1 and >0 we can conclude that the orbit is a retrogade elliptical orbit\n',e_0,magne_0)
fprintf('The obatined Orbital distance from the Earth center at Perigee is %0.2fkm\n',r_P)
fprintf('The obatined Orbital distance from the Earth center at Apogee is %0.2fkm\n',r_A)
fprintf('The obatined semi major axis is %0.2fkm\n',a)
fprintf('The obatined Orbital Period is %0.2fs\n',r_P)


%B)

% Simulation (ODE Solver)
Options = odeset('MaxStep', 10, 'RelTol', 1e-11, 'AbsTol', 1e-11);
[Time, States] = ode45(@TBP, t_initial:t_step:t_final, [r_0; v_0], Options);

% Output
% Position Vector [km]
r_x = States(:,1);
r_y = States(:,2);
r_z = States(:,3);
% Velocity Vector [km]
v_x = States(:,4);
v_y = States(:,5);
v_z = States(:,6);

figure (1)

plot3(r_x, r_y, r_z,"m", 'LineWidth', 2);
hold on
plot3(r_x(1), r_y(1), r_z(1), 'r*', 'LineWidth', 2);
hold on
% Adding the Earth

[X,Y,Z] = sphere;
Earth = surf(Re*X, Re*Y, -Re*Z);
alpha = 1; % Transparency of the Surface of the Earth
cdata = imread('Earth_Texture.jpg');
set(Earth, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
hold on

title('Spacecraft''s Trajectory in the Earth-Centered Inertial Frame','Interpreter','Latex');
xlabel('$\hat{x}$ [km]','Interpreter','Latex');
ylabel('$\hat{y}$ [km]','Interpreter','Latex');
zlabel('$\hat{z}$ [km]','Interpreter','Latex');
legend('Spacecraft''s Trajectory','Initial Condition','Interpreter','Latex')

axis equal
grid on-    
grid minor
set(gcf,'Color','W');
fontsize(gcf, 24, 'Points')
fontname(gcf,'Times New Roman')


%C) a.
theta = (0:1:360) ;
r_parametric_solution = ((magnh_0^2)/mui)./(1 + (magne_0*cosd(theta))) ; %make sure to use the . before the PEMDAS operations 
r_parametric_solution_vector  = (r_parametric_solution.*(cosd(theta)).*(e_0/magne_0) )+(r_parametric_solution.*sind(theta).*(cross((h_0/magnh_0),(e_0/magne_0)))) ; %Perifocal coordinate system


plot3(r_parametric_solution_vector(1,:),r_parametric_solution_vector(2,:),r_parametric_solution_vector(3,:),'k','LineWidth',2)
hold on 

%D) a.
v_circ_perigee = ((mui/r_P)^(1/2)) ; %computing the velocity of the spacecraft at perigee to be in a circular orbit
fprintf('The velocity of the spacecraft at perigee in a circular orbit is %0.4fkm/s\n',v_circ_perigee)

%a_p_
T_P = (2*pi)*(mui^(-1/2))*(r_P^(3/2)) ;

r_P_vector = (r_P.*(e_0/magne_0)) ;
v_circ_perigee_vector =  (v_circ_perigee.*(cross((e_0/magne_0),(h_0/magnh_0)))) ;

Options = odeset('MaxStep', 10, 'RelTol', 1e-11, 'AbsTol', 1e-11);
[Timep, Statep] = ode45(@TBP, t_initial:t_step:T_P, [r_P_vector;v_circ_perigee_vector], Options);

plot3(Statep(:,1), Statep(:,2),Statep(:,3),'b','LineWidth', 2);
hold on

%E) a.
v_circ_apogee = ((mui/r_A)^(1/2)) ;  %computing the velocity of the spacecraft at apogee to be in a circular orbit
fprintf('The velocity of the spacecraft at apogee in a circular orbit is %0.4fkm/s\n',v_circ_apogee)

r_A_vector = (r_A.*(e_0/magne_0)) ;
v_circ_apogee_vector =  (v_circ_apogee.*(cross((e_0/magne_0),(h_0/magnh_0)))) ;
T_A = (2*pi)*(mui^(-1/2))*(r_A^(3/2));
Options = odeset('MaxStep', 10, 'RelTol', 1e-11, 'AbsTol', 1e-11);
[Timea, Statea] = ode45(@TBP, t_initial:t_step:T_A, [r_A_vector;v_circ_apogee_vector], Options);

plot3(Statea(:,1),Statea(:,2), Statea(:,3),'g','LineWidth', 2);
hold on





%TBP FUNCTION

function State_dot = TBP(~,State)

% Constants of the Problem
G = 6.67259e-20; % Gravitational Constant [km^3/(kg*s^2)] 
Me = 5.9722e24; % Mass of the Earth [kg]
muI = G*Me; % Gravitational Parameter of the Earth [km^3/s^2]
Re = 6371; % Mean Radius of the Earth [km] 
J2 = 0.00108263; % Second Zonal Harmonic of the Earth [Unitless]

% Dynamics

% Position Vector [km]
x = State(1);
y = State(2);
z = State(3);

r_vector = [x; y; z];

r = norm(r_vector);

% Velocity Vector [km/s]
x_dot = State(4);
y_dot = State(5);
z_dot = State(6);

v_vector = [x_dot; y_dot; z_dot];
% Dynamics
r_dot_vector = v_vector;
v_dot_vector = -muI/r^3*r_vector;

% Output
State_dot = [r_dot_vector; v_dot_vector];

end