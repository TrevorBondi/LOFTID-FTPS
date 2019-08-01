function[L, A, r_in, r_out, d_frac, angle_c_rad, angle_s_rad] = FTPS_Mass(n, d1, d2, dp, ply, angle_c_deg, angle_s_deg_vertical, c_factor, nose_geometry, SiC_density)
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
angle_c_rad = deg2rad(angle_c_deg);
angle_s_rad_vertical = deg2rad(angle_s_deg_vertical);
angle_s_rad = -pi/2+angle_c_rad+angle_s_rad_vertical;
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

%fprintf('Function Call Success\n')
end