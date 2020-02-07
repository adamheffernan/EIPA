%Adam Heffernan 100977570 2/7/2020

q_0 = 1.60217653e-19;             % electron charge
hb = 1.054571596e-34;             % Dirac constant
h = C.hb * 2 * pi;                    % Planck constant
m_0 = 9.10938215e-31;             % electron mass
kb = 1.3806504e-23;               % Boltzmann constant
eps_0 = 8.854187817e-12;          % vacuum permittivity
mu_0 = 1.2566370614e-6;           % vacuum permeability
c = 299792458; % speed of light

nx = 50;
ny = 50;

G = sparse(nx*ny, ny*nx);
alpha = (C.hb^2) / (2 * C.m_0);

mapping_equation = @(i,j) j + (i - 1)*ny;


for i=1:nx
    for j=1:ny
        current_cell = mapping_equation(i,j);
        nx_left = mapping_equation(i-1,j);
        nx_right = mapping_equation(i+1,j);
        ny_bottom = mapping_equation(i,j-1);
        ny_top = mapping_equation(i,j+1);
        
        if (i == 1 || i == nx)
            G(current_cell,current_cell) = 1; 
        elseif (j == 1 || j == ny)
            G(current_cell,current_cell) = 1;
            
        elseif (i > 10 && i < 20 && j > 10 && j < 20)
            
            G(current_cell,current_cell) = -2;
        else
            G(current_cell,current_cell) = -4;
            G(current_cell,nx_left) = 1;
            G(current_cell,nx_right) = 1; 
            G(current_cell,ny_bottom) = 1; 
            G(current_cell,ny_top) = 1;
        end
    end
end
figure(2);
spy(G);
E_plot = zeros(nx,ny,9);

figure(3);
[E,D] = eigs(G,9,'SM');
plot(D);
figure(1)
hold on;

for current_cell=1:9
    for i=1:nx
        for j=1:ny
            E_plot(i,j,current_cell) = E(mapping_equation(i,j),current_cell);
        end
    end
    subplot(3,3,current_cell);
    surf(E_plot(:,:,current_cell));
end


