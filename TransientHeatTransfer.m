%% Heat Transfer of Turbulent Boundary Layer in Accelerating Combustion Gases
%
% Heat Transfer Fundamental Eq:
% T1(t+1) = (2 * h * dt)/(rho * cp * dx) * (Tg - T1(t)) +
% (2*k/rho/cp*dt/dx^2)*(T2(t) - T1(t)) + T1(t)
%


% inputs
thickness = 0.2; % [inches]
hg = 2584.8; % [W/m2/degC] heat transfer coeff of combustion gases
ha = 50; % [W/m2/degC] heat transfer coeff of ambient air
k = 0.6918; % [W/m/degC]
T_start = 298; % starting temperature (

NoFiniteElems = 50; % number of finite elements
dt = 0.0001; % time increment [s]
dT = 10; % total time simulated [s]

rho = 8000; % [kg/m3] density of material
cp = 490; % [J/kg/degK] heat capacity of chamber material. Ref: Cp_steel = 490
Tg = 2448; % [degK] combustion product flame temp



%% 1D transiet heat transfer
% unit conversions
t = thickness * 0.0254;

dx = t / NoFiniteElems;
T = ones(NoFiniteElems, dT/dt) .* T_start;

A = zeros(NoFiniteElems, NoFiniteElems);
A(1,1) = -2*hg*dt/rho/cp/dx - 2*k*dt/rho/cp/dx^2 + 1;
A(1,2) = 2*k*dt/rho/cp/dx^2;

A(end,end) = 2*dt/rho/cp/dx*(-k/dx - ha) + 1;
A(end,end-1) = 2*k*dt/rho/cp/dx^2;

B = zeros(NoFiniteElems,1);
B(1,1) = 2*hg*dt/rho/cp/dx*Tg;
B(end,1) = 2*ha*dt*T_start/rho/cp/dx;

for i = 2:NoFiniteElems-1
    A(i,i) = -2*k*dt/rho/cp/dx^2 + 1;
    A(i,i-1) = k*dt/rho/cp/dx^2;
    A(i,i+1) = k*dt/rho/cp/dx^2;
end

for i = 2:(dT/dt)
    T(:,i) = A * T(:,i-1) + B;
end


%% Plotting
for ind = [2000 5000 10000 20000 50000 100000]
    Tsec = T(:,ind);
    time = ind/length(T) * dT;
    y = 0:(thickness/NoFiniteElems):thickness;
    x = 0:1/NoFiniteElems:1;
    
    [X, Y]=meshgrid(x,y);
    Z = ones(length(x),length(y));
    for i = 1:NoFiniteElems
        Z(i,:) = Tsec(i);
    end
    figure()
    surf(X,Y,Z);
    
    view(0,90)
    title("T = " + time + "s")
    xlabel("x")
    ylabel("y distance from wall [in]")
    
    c = colorbar;
    c.Label.String = 'Temperature [K]';
%     caxis([T_start,Tsec(end)]);
end

%% 
[X,Y] = meshgrid(1:0.5:10,1:20);
Z = sin(X) + cos(Y);
surf(X,Y,Z)
view(0,90)