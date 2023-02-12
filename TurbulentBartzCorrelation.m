%% equilibrium reaction properties
A_inch2 = 3.10368; % chosen area
M = 0.207; % mach number
m = 1; % temperature dependence of viscosity

Astar_inch2 = 0.38765943; % throat area
Dstar_inch = 0.70282665; % throat diameter
mu_milipoise = 0.85653; % dynamic viscosity
Pr = 0.5004; % prandtl number
Cp_kJ_kg = 2.1795; % spec heat @ const. pressure
Po_psi = 400; % chamber pressure
g = 9.8; % gravity constant
Cstar = 1531.7; % characteristic velocity m/s
rc_inch = 0.4; % wall radius of curvature @ throat
gamma = 1.2201; % specific heat ratio
Tw = 298; % starting wall temp degK
To = 298; % stagnation temp degK

sigma = 1 / (0.5 * Tw / To * (1 + (gamma - 1)/2*M^2) + 1/2)^(0.8 - m/5)/(1 + (gamma - 1)/2*M^2)^(m/5);

% throat params (not used)
A_inch2_throat = Astar_inch2;
gamma_throat = 1.2346;
pressure_throat = 400 / 1.7936; % psi
M_throat = 1;
Cp_kJ_kg_throat = 1.9904;
mu_milipoise_throat = 0.79525;
Pr_throat = 0.5488;





% unit conversions
A = A_inch2 * 0.00064516;
Astar = Astar_inch2 * 0.00064516;
Dstar = Dstar_inch * 0.0254;
mu = mu_milipoise * 0.0001;
Cp = Cp_kJ_kg * 1000;
Po = Po_psi*6894.757293;
rc = rc_inch * 0.0254;

A_inch2_throat = A_inch2_throat*0.00064516;
mu_milipoise_throat = mu_milipoise_throat * 0.0001;
Cp_kJ_kg_throat = Cp_kJ_kg_throat * 1000;
pressure_throat = pressure_throat * 6894.757293;

%% get yvals from contour
avals = yvals .^2 *pi;
avals = avals .*0.00064516;
hg = 0.026 ./ Dstar.^0.2 .* (mu.^0.2 .* Cp./Pr.^0.6) .* (Po.*g./Cstar).^0.8 .* (Dstar./rc).^0.1 .* (Astar./avals).^0.9 .* sigma;

figure()
plot(xvals./xvals(end), hg)
title('Heat Transfer Coeff over Engine')
ylabel("h_g [W/m2/degC]")
xlabel('X/L')


% figure()
% contur = [xvals./xvals(end);hg]';
% contourf(contur)

figure()
[X, Y] = meshgrid(xvals./xvals(end),yvals);
HG = [];
for i = 1:length(xvals)
    HG(:,i) = hg(i).*ones(length(xvals),1);
end

z = [HG];
[c,h] = contourf(z);
colorbar
legend('W/m2/degC')

%%

