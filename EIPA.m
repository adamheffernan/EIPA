global C;
C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458; % speed of light

nx = 75;
ny = 50;

dx = 1;
dy = 1;
G = sparse(nx*ny, ny*nx);
alpha = (C.hb^2) / (2 * C.m_0);

map = @(i,j) j + (i - 1)*ny;


for i=1:nx
    for j=1:ny
        n = map(i,j);
        nxm = map(i-1,j);
        nxp = map(i+1,j);
        nym = map(i,j-1);
        nyp = map(i,j+1);
        
        if (i == 1 || i == nx)
            G(n,n) = 1 / dx^2;
        elseif (j == 1 || j == ny)
            G(n,n) = 1 / dy^2;
            
        elseif (i > 10 & i <20 & j > 10 & j <20)
            G(n,n) = -1 / dx^2 + -1 / dy^2;
            G(n,nxm) = 1 / dx^2;
            G(n,nxp) = 1 / dx^2;
            G(n,nym) = 1 / dy^2;
            G(n,nyp) = 1 / dy^2;
        else
            G(n,n) = -2 / dx^2 + -2 / dy^2;
            G(n,nxm) = 1 / dx^2;
            G(n,nxp) = 1 / dx^2;
            G(n,nym) = 1 / dy^2;
            G(n,nyp) = 1 / dy^2;
        end
    end
end

%spy(G);
[E,D] = eigs(G,9,'SM');

E_plot = zeros(nx,ny,9);
for n=1:9
    for i=1:nx
        for j=1:ny
            E_plot(i,j,n) = E(map(i,j),n);
        end
    end
end

f = figure();
hold on;
for n=1:9
    subplot(3,3,n);
    surf(E_plot(:,:,n),'EdgeColor','none','LineStyle','none','FaceLighting','phong');
end