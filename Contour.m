% %% conical 
% chamber_length =
% converge_angle

%% parabolic
chamber_length = 5.68182; % inch
chamber_diameter = 1.98789; % inch
throat_diameter = 0.70283; % inch
exit_diameter = 1.43016; % inch
e = (exit_diameter / throat_diameter)^2;
nozzle_length = 0.8 * ((sqrt(e) - 1)*throat_diameter/2)/tand(15);
theta_i = 45; % degrees

rc1 = 1.5 *  throat_diameter/2; % radius of curvature at throat inlet
rc2 = 0.4 * throat_diameter / 2; % radius of curvature at throat exit

syms x
circle1 = sqrt(rc1^2 - (x-(chamber_length))^2)+chamber_diameter/2-rc1; % ( x - h )^2 + ( y - k ) ^2  = r ^2
xl = chamber_length + rc1*cosd(45);
yl = chamber_diameter/2 -rc1 + rc1*cosd(45);
drop = abs((chamber_diameter/2 - throat_diameter/2) - 2*(rc1 - rc1*cosd(45)));
b = yl + xl;
line = -x + b;
circle2 = -sqrt(rc1^2 - (x-(chamber_length + 2*rc1*cosd(45) + drop))^2) + chamber_diameter/2-rc1+rc1*cosd(45) - drop +rc1*cosd(45);
circle3 = -sqrt(rc2^2 - (x - (chamber_length + 2*rc1*cosd(45) + drop))^2) + chamber_diameter/2 - rc1+2*rc1*cosd(45) - drop - rc1 + rc2;

y1 = chamber_diameter/2 - rc1 +2*rc1*cosd(45)- rc1 + rc2 - rc2*cosd(30) - drop; % start of parabola
x1 = chamber_length + 2*rc1*cosd(45) + rc2*cosd(60) + drop; % start of parabola
y2 = exit_diameter/2;
x2 = chamber_length + 2*rc1*cosd(45) + drop + nozzle_length;

a = (y2-y1+tand(30)*(x1-x2))/(x2^2 - 2*x1*x2 + x1^2);
c = y1 + a*x1^2 - tand(30)*x1;
b = tand(30) - 2*a*x1;
syms parabola(x)
parabola(x) = a*x^2 + b*x + c; 
syms y(x)
y(x) = piecewise((x < chamber_length) & (x >= 0), chamber_diameter/2, (x >= chamber_length) & (x < chamber_length + rc1*cosd(45)), circle1, (x >= chamber_length + rc1*cosd(45)) & (x < chamber_length + rc1*cosd(45) + drop), line, (x >= chamber_length + rc1*cosd(45) + drop) & (x < chamber_length + 2*rc1*cosd(45) + drop), circle2, (x >= chamber_length + 2*rc1*cosd(45) + drop) & (x < chamber_length + 2*rc1*cosd(45) + drop + rc2*cosd(60)), circle3, (x >= x1) & (x < x2), parabola);
xvals = 0:0.01:x2;
yvals = double(y(xvals));
%% plotting
figure()
plot(xvals,yvals)
axis equal
ylim([0,2])
xlim([0,x2+1])
