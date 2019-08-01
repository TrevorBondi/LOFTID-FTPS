% NASA LOFTID F-TPS Mass Modeler Prototype
clc; clear all; close all; hold off;
fprintf('Diameters of a given unit yields areas of unit^2\n')

%% HIAD Characteristics
ply = 4;
n = 7;
d1 = 15.7;
d2 = 7.697;
dp = 49;
angle_s_rad_vertical = deg2rad(15);
angle_c_rad = deg2rad(70);
angle_s_rad = -pi/2+angle_c_rad+angle_s_rad_vertical;
nose_geometry = 's';
SiC_density = 12.7; %oz/yd^2
c_factor = 1.20;

%{
KEY
n = number of tori, including smaller shoulder torus
d1 = diameter of normal tori
d2 = diameter of shoulder torus
dp = diameter of enclosure/payload structure
angle_s_rad_vertical = angle F-TPS extends over the last torus relative to vert.
angle_s_rad = angle the F-TPS extends over the last torus relative to
line normal to half cone angle.
angle_c_rad = half cone angle of the F-TPS skirt (w/ respect to vertical)
nose_geometry = 's' (sphere) or 'c' (cone) - other strings return no geometry
c_factor = correction factor accounts for seams and other manufacturing factors
%}

%% Pre-calculations to save time/space
d_frac = (d1-d2)/(d1+d2);
r_in = zeros(1,6);
r_out = zeros(1,6);

%% Length of F-TPS half cross section
L(1) = (n-1)*d1+0.5*(d1+d2)*cos(asin(d_frac))-d1/2;
L(2) = d2/2*(pi-2*asin(d_frac));
L(3) = 0.5*(d1+d2)*cos(asin(d_frac));
L(4) = d1/2*(angle_s_rad+2*asin(d_frac));
L(5) = d1/2*(1+tan(pi/4-angle_c_rad/2))/tan(angle_c_rad);

%{
"L" is a vector containing 5 values representing the lengths that make up the
cross section of the F-TPS.

The lengths L(1) through L(5) are named non sequentially. That is to say, L(1)
represents the length running beneath the tori that connects with the outer
shoulder tangentially. L(5) is the length of the gap between L(1) and the
body/wall of the payload or enclosure. L(2) is the arc length of the shoulder,
L(4) is the arc length of the inner most shoulder that wraps slightly around the
last large torus, and L(3) is the flat spot created between these two shoulders.
This order arose because of the order the geometry was derived.
%}

%% Surface Area of cones and shoulders using frustum approximation (A = 2 * pi * 0.5 * (outer radius + inner radius) * length of curve)
if nose_geometry == 's'
    A(1) = pi*((dp/2)^2+(dp/2*(1/cos(angle_c_rad)-tan(angle_c_rad)))^2); % Exact spherical nose area
elseif nose_geometry == 'c'
    A(1) = pi*(dp/2)^2/sin(angle_c_rad); % Exact conical nose area
else
    A(1) = 0;
end

A(2) = pi*(dp+(L(1)+L(5))*sin(angle_c_rad))*(L(1)+L(5)); % Exact front facing cone area

r_out(3) = dp/2+(L(1)+L(5))*sin(angle_c_rad)-d2/2*sin(pi/2-angle_c_rad)+d2/2*sin(2*asin(d_frac)-pi/2+angle_c_rad);
r_in(3) = dp/2+(L(1)+L(5))*sin(angle_c_rad);
A(3) = pi*(r_in(3)+r_out(3))*L(2); % Approximate outer shoulder surface area (approximation)

r_out(4) = dp/2+(L(1)+L(5))*sin(angle_c_rad)-d2/2*sin(pi/2-angle_c_rad)+d2/2*sin(2*asin(d_frac)-pi/2+angle_c_rad);
r_in(4) = r_out(4)-L(3)*cos(2*asin(d_frac)-pi/2+angle_c_rad);
A(4) = pi*(r_out(4)+r_in(4))*L(3); % Exact surface area of flat spot between shoulders

r_out(5) = dp/2+(L(1)+L(5))*sin(angle_c_rad)-d2/2*sin(pi/2-angle_c_rad)+d2/2*sin(2*asin(d_frac)-pi/2+angle_c_rad)-L(3)*cos(2*asin(d_frac)-pi/2+angle_c_rad);
r_in(5) = dp/2+(L(5)+(n-1.5)*d1)*sin(angle_c_rad)-d1/2*cos(angle_c_rad) - d1/2*sin(pi/2+angle_s_rad-angle_c_rad);
A(5) = pi*(r_in(5) + r_out(5))*L(4); % Approximate inner shoulder surface area

r_in(6) = dp/2+(L(1)+L(5))*sin(angle_c_rad);
r_out(6) = r_in(6)+d2/2*(1-cos(angle_c_rad));
A(6) = pi*(r_in(6)+r_out(6))*d2/2*angle_c_rad;

% fprintf('Nose Area: %f\n', A(1))
% fprintf('Front Conical Area: %f\n', A(2))
% fprintf('Outer Shoulder: %f\n', A(3))
% fprintf('Rear Conical Area: %f\n', A(4))
% fprintf('Inner Shoulder: %f\n', A(5))
% fprintf('Outer Ply Area: %f\n', sum(A(1:5)))
% fprintf('Inner Ply Area %f\n', A(1)+A(2)+A(6))
A_outer = sum(A(1:5));
A_inner = A(1)+A(2)+A(6);
A_total = A_outer+(ply-1)*A_inner;
A_corrected = c_factor*A_total;
SiC_mass = A_corrected/36^2*SiC_density;

%{
"A" is a vector containing 5 values representing the surface areas of each
geometric component of the F-TPS (nose, front-facing cone, outer shoulder,
inner shoulder, flat spot between shoulders).

The areas produced by the above calculations are sequentially ordered.
That is to say, the first element of A corresponds to the nose surface area,
the second element corresponds to the forward facing conical surface of the
F-TPS, the third element corresponds to the area of the outer shoulder,
the fifth element corresponds to the inner shoulder, and the fourth element
corresponds to the flat spot (rear conical area) created between these
two shoulders.
%}

%% Plotting Cross Section as Visual Aid
figure(1); hold on; set(gca,'DefaultLineLineWidth',3); axis equal; grid on; title('LOFTID HIAD F-TPS Cross Section')

P1 = [dp/2,0];
P2 = [dp/2+(L(1)+L(5))*sin(angle_c_rad),(L(1)+L(5))*cos(angle_c_rad)];
plot([P1(1), P2(1)], [P1(2), P2(2)], 'r') % Plots front conical cross section
plot(-[P1(1), P2(1)], [P1(2), P2(2)], 'r')

y_range = [P2(2):0.1:(L(1)+L(5))*cos(angle_c_rad)+d2/2*sin(angle_c_rad)+d2/2*cos(2*asin(d_frac)-pi/2+angle_c_rad)];
x_offset = dp/2+(L(1)+L(5))*sin(angle_c_rad)-d2/2*cos(angle_c_rad);
y_offset = (L(1)+L(5))*cos(angle_c_rad)+d2/2*sin(angle_c_rad);
x_function = sqrt((d2/2)^2-(y_range-y_offset).^2)+x_offset;
x_reflect = -x_function;
plot(x_function, y_range, 'r') % Plots outer shoulder
plot(x_reflect, y_range, 'r')

x_offset = dp/2+(L(5)+(n-1.5)*d1)*sin(angle_c_rad)-d1/2*cos(angle_c_rad);
y_offset =(L(5)+(n-1.5)*d1)*cos(angle_c_rad)+d1/2*sin(angle_c_rad);
x_range = [x_offset-d1/2*sin(pi/2-angle_c_rad+angle_s_rad):0.1:r_in(4)];
x_reflect = -x_range;
y_function = sqrt((d1/2)^2-(x_range-x_offset).^2)+y_offset;
plot(x_range, y_function, 'r') % Plots inner shoulder
plot(x_reflect, y_function, 'r')

P3 = [r_out(4), max(y_range)];
P4 = [r_in(4), y_function(end)];
plot([P3(1), P4(1)], [P3(2), P4(2)], 'r') % Plots flat between shoulders
plot(-[P3(1), P4(1)], [P3(2), P4(2)], 'r')

if nose_geometry == 's'
    x_range = [0:0.1:dp/2];
    x_reflect = -x_range;
    y_offset = dp/2*tan(angle_c_rad);
    r = dp/(2*cos(angle_c_rad));
    y_function = -sqrt(r^2-x_range.^2)+y_offset;
    plot(x_range, y_function, 'b')
    plot(x_reflect, y_function, 'b')
elseif nose_geometry == 'c'
    P5 = [0,-dp/2*tan(pi/2-angle_c_rad)];
    plot([P5(1),P1(1)],[P5(2),P1(2)], 'b')
    plot(-[P5(1),P1(1)],[P5(2),P1(2)], 'b')
end
%% Cosmetic Plot of Torus cross sections


for i = 1:(n-1)
    hold on;
    x_cent = dp/2+(L(5)+(i-0.5)*d1)*sin(angle_c_rad)-d1/2*cos(angle_c_rad);
    y_cent = (L(5)+(i-0.5)*d1)*cos(angle_c_rad)+d1/2*sin(angle_c_rad);
    r = d1/2;
    circle(x_cent, y_cent, r);
    circle(-x_cent, y_cent, r);
end

x_cent = dp/2+(L(1)+L(5))*sin(angle_c_rad)-d2/2*cos(angle_c_rad);
y_cent = (L(1)+L(5))*cos(angle_c_rad)+d2/2*sin(angle_c_rad);
r = d2/2;
circle(x_cent, y_cent, r);
circle(-x_cent, y_cent, r);

%% Output to Command Window
fprintf('Outer Ply Surface Area: %0.1f yd^2\n', A_outer/36^2);
fprintf('Inner Ply Surface Area: %0.1f yd^2\n', A_inner/36^2);
fprintf('Uncorrected Total Area: %0.1f yd^2\n', A_total/36^2);
fprintf('Corrected Total Area: %0.1f yd^2\n', A_corrected/36^2);
fprintf('SiC Fabric Mass: %0.1f oz\n', SiC_mass);

function circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.1:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp,'--k');
end