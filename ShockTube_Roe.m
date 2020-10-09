%% 1D Shock Tube using Roe and Adjustable Time Stepping
close all
clear
clc
%% Initialization of Parameters
N       = 800 ;             % Number of grid points
gamma   = 1.4 ;
endTime = 0.2 ;             % Physical End Time
CFL     = 1.0 ;
%% Grid Initialization
x       = linspace(0,1,N+1) ;
xc      = 0.5 * ( x(1:N) + x(2:N+1) ) ;
xc(1)   = 0 ;               % Xmin
xc(N)   = 1 ;               % Xmax
time    = 0 ;
%% Initial Conidtions
denistyR    = 0.125 ;       densityL    = 1.0 ;
pressureR   = 0.1   ;       pressureL   = 1.0 ;

rho     = zeros(N,1) ;
p       = zeros(N,1) ;
u       = zeros(N,1) ;
Flux    = zeros(N,3) ;

Flux(1,2) = 1 ;         % momentum left bc
Flux(N,2) = 0.1 ;       % momentum right bc

for i =  1:N
    if i<=N/2
        rho(i)  = densityL  ;
        p(i)    = pressureL ;
    else
        rho(i)  = denistyR  ;
        p(i)    = pressureR ;
    end
end
e   = p/(gamma-1) + 0.5*rho.*u.*u ;
h   = gamma/(gamma-1)*p./rho + 0.5*u.*u ;
%%
new_rho = rho ;
new_u   = u   ;
new_e   = e   ;

while time <= endTime
    
    for i=2:N-1
        p(i)    = (gamma-1)*(e(i) - 0.5*rho(i)*u(i)*u(i)) ;
        h(i)    = gamma/(gamma-1)*p(i)/rho(i) + 0.5*u(i)*u(i) ;
    end
    a       = sqrt(gamma*p./rho) ;
    lamda   = max(a) ;
    max_vel = max(u) ;
    
    dt      = CFL/N/(max_vel+lamda) ;  % adjustable Time Step
    time    = time + dt ;
    
    for i=2:N-1
        dx          = xc(i) - xc(i-1) ;
        
        mom_R       = rho(i+1)*u(i+1) ;     rho_R = rho(i+1) ;      u_R = u(i+1) ;      p_R = p(i+1) ;
        mom_L       = rho(i-1)*u(i-1) ;    	rho_L = rho(i-1) ;      u_L = u(i-1) ;      p_L = p(i-1) ;
        
        e_R = e(i+1) ;      h_R = h(i+1) ;
        e_L = e(i-1) ;      h_L = h(i-1) ;
        
        rho_RL  = sqrt(rho_R*rho_L) ;
        u_RL    = ( sqrt(rho_R)*u_R + sqrt(rho_L)*u_L )/( sqrt(rho_R) + sqrt(rho_L) ) ;
        h_RL    = ( sqrt(rho_R)*h_R + sqrt(rho_L)*h_L )/( sqrt(rho_R) + sqrt(rho_L) ) ;
        a_RL    = sqrt( (gamma-1)*(h_RL - 0.5*u_RL*u_RL) ) ;
        
        lamda1  = u_RL ;
        lamda2  = u_RL + a_RL ;
        lamda3  = u_RL - a_RL ;
        
        drho    = rho_R - rho_L ;
        du      = u_R   - u_L   ;
        dp      = p_R   - p_L   ;
        
        dv1     = drho - dp/a_RL/a_RL ;
        dv2     = du + dp/rho_RL/a_RL ;
        dv3     = du - dp/rho_RL/a_RL ;
        
        vel_flux_R  = rho_R*u_R*u_R +p_R ;
        vel_flux_L  = rho_L*u_L*u_L +p_L ;
        
        energy_flux_R   = u_R * ( e_R + p_R ) ;
        energy_flux_L   = u_L * ( e_L + p_L ) ;
        
        rho_flux     = mom_R - max(0,lamda1)*dv1 - 0.5*rho_RL/a_RL*( max(0,lamda2)*dv2 - max(0,lamda3)*dv3 )  ;
        
        vel_flux     = vel_flux_R - u_RL*max(0,lamda1)*dv1 - 0.5*rho_RL/a_RL*( max(0,lamda2)*lamda2*dv2 - max(0,lamda3)*lamda3*dv3 )  ;
        
        energy_flux  = energy_flux_R - u_RL*u_RL*max(0,lamda1)*dv1 - 0.5*rho_RL/a_RL*( max(0,lamda2)*(h_RL+a_RL*u_RL)*dv2 - max(0,lamda3)*(h_RL-a_RL*u_RL)*dv3 )  ;
        
        Flux(i,1) = rho_flux ;
        Flux(i,2) = vel_flux ;
        Flux(i,3) = energy_flux ;
        
    end
    
    for i=2:N-1
        
        new_rho(i)  = rho(i)        - 0.5*dt/dx * ( Flux(i+1,1) - Flux(i-1,1) ) ;
        vel_flux    = rho(i)*u(i)   - 0.5*dt/dx * ( Flux(i+1,2) - Flux(i-1,2) ) ;
        new_e(i)    = e(i)          - 0.5*dt/dx * ( Flux(i+1,3) - Flux(i-1,3) ) ;
        
        new_u(i)    = vel_flux/new_rho(i) ;
        
    end
    
    rho     = new_rho ;
    u       = new_u ;
    e       = new_e ;
    
end

pressure = dlmread('pressure.dat') ;
density  = dlmread('density.dat')  ;
velocity = dlmread('velocity.dat') ;

figure(1)
hold on
plot(density(:,1),density(:,2),'r--','MarkerSize',12,'MarkerIndices',1:10:length(velocity),'LineWidth',2);
plot(xc,rho,'k','LineWidth',2);
xlabel(' x ','FontSize',18,'FontWeight','bold');
ylabel(' Density ','FontSize',18,'FontWeight','bold');
legend('Analytical','Roe','Location','northeast','FontSize',16,'FontWeight','bold');
%print(gcf,'Density_Roe.jpg','-dpng','-r300');

figure(2)
hold on
plot(pressure(:,1),pressure(:,2),'r--','MarkerSize',12,'MarkerIndices',1:10:length(velocity),'LineWidth',2);
plot(xc,p,'k','MarkerSize',12,'MarkerIndices',1:10:141,'LineWidth',2);
xlabel(' x ','FontSize',18,'FontWeight','bold');
ylabel(' Pressure ','FontSize',18,'FontWeight','bold');
legend('Analytical','Roe','Location','northeast','FontSize',16,'FontWeight','bold');
%print(gcf,'Pressure_Roe.jpg','-dpng','-r300');

figure(3)
hold on
plot(velocity(:,1),velocity(:,2),'r--','MarkerSize',12,'MarkerIndices',1:10:length(velocity),'LineWidth',2);
plot(xc,u,'k','MarkerSize',12,'MarkerIndices',1:10:141,'LineWidth',2);
xlabel(' x ','FontSize',18,'FontWeight','bold');
ylabel(' Velocity ','FontSize',18,'FontWeight','bold');
legend('Analytical','Roe','Location','south','FontSize',16,'FontWeight','bold');
%print(gcf,'Velocity_Roe.jpg','-dpng','-r300');