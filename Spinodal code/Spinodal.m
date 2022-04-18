%this is a code to simulate ferromagnetic domain growth
%we are going to decide the initial condititon
% all the simulation parameters are non-dimensional

function[]=PF_mag(n,tot_tsteps)

%specifying the input parameters 
%tot_tsteps = the total no. of timesteps the simulation has to run  
%n = the total no. of grid points along the square edges of a system
%dx = the grid spacing
%dt = timestep size

%declaring a few variables as global
global dx
global dy
global H_b
global kappa
global sqdim

sqdim=n; %the grid dimensions are stored here  
dx=1; %the grid spacing along x
dy=1; %the grid spacing along y
dt=0.01; %the timestep size
H_b=1.0; %the parameter which sets the height of the activation barrier
kappa=1.0; %the gradient energy density coeffcient
M=1.0; %the mobility of the solutes 

%creating the initial condition 
c = 0.5 + 0.05*(2.0*randn(n) - 1);
          
%storing the initial profile 
c_init = c;

%solving the Cahn-Hilliard equation
%starting the time loop
for t=1:tot_tsteps
    %the space loops
    for i=1:n 
        for j=1:n
            %calculating the laplacian of mu
            mu_current = mu_calc(c,i,j);
            mu_right   = mu_calc(c,pbc(i+1,n),j); 
            mu_left    = mu_calc(c,pbc(i-1,n),j);
            mu_up      = mu_calc(c,i,pbc(j+1,n));
            mu_down    = mu_calc(c,i,pbc(j-1,n));
            
            lapl_mu = ((mu_left-(2.0*mu_current)+mu_right)/(dx^2)) + ((mu_up - (2.0*mu_current) + mu_down)/(dy^2));
            
            %calculating the update value of c
            update_c(i,j) = M*dt*lapl_mu; 
                   
        end
    end   
 
    %restarting the spacee loops to update the c field
    for i=1:n
        for j=1:n
            c(i,j) = c(i,j) + update_c(i,j);
        end
    end     
end    
    
%viewing the order parameter profile
plotres(c_init,c);

%imagesc(phi);
%colormap(gray)
%axis equal
%axis off




