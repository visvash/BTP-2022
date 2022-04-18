%this is a code to simulate ferromagnetic domain growth
%we are going to decide the initial condititon
% all the simulation parameters are non-dimensional

function[]=PF_mag(c,tot_tsteps,n,h)

%specifying the input parameters 
%c = the choice of initial conditions; c=1 planar interface, c=2, circular
%domain
%tot_tsteps = the total no. of timesteps the simulation has to run  
%n = the total no. of grid points along the square edges of a system
%dx = the grid spacing
%dt = timestep size
%h = the externally imposed magnetic field
 
dom_rad=40; %initializing the domain radius in the event of starting with a circular domain 
dx=0.5; %the grid spacing along x
dy=0.5; %the grid spacing along y
dt=0.01; %the timestep size
H_b=1.0; %the parameter which sets the height of the activation barrier
kappa=1.0; %the gradient energy density coeffcient
Gamma=1.0; %the mobility of the interface 

%creating the initial condition 
if c==1 %creating a planar interface
    for i=1:n
        for j=1:n
           if j >= (n/2)
               phi(i,j) = 1.0;
           else 
               phi(i,j) = -1.0;
           end
        end
    end 
elseif c==2 %creating a circular domain    
        for i=1:n
            for j=1:n
                if  (sqrt((i-(n/2))^2 + (j -(n/2))^2)) < dom_rad
                    phi(i,j) = -1.0; 
                else 
                    phi(i,j) = 1.0;
                end
            end 
        end
end   

%storing the initial profile 
phi_init = phi;
phi_init_profile = phi(n/2,:);

%solving the Allen-Cahn equation
%starting the time loop
for t=1:tot_tsteps
    %the space loops
    for i=1:n 
        for j=1:n
            %evaluating the derivative of bulk energy density
            dfb_dphi = -4.0*H_b*phi(i,j)*(1.0-phi(i,j)^2);
            %evaluating the term due to the gradient free energy
            lapl_phi = ((phi(pbc(i+1,n),j) - (2.0*phi(i,j)) + phi(pbc(i-1,n),j))/(dx^2)) + ((phi(i,pbc(j+1,n)) - (2.0*phi(i,j)) + phi(i,pbc(j-1,n)))/(dy^2));
            term_grad_free = kappa*lapl_phi;
            %evaluting the driving force term
            dg_dphi = (3.0*(1.0-phi(i,j)^2))/2.0;
            driv_force= h*dg_dphi;
            %calculating the value by which phi would be updated 
            update_phi(i,j) = -Gamma*dt*(dfb_dphi - term_grad_free - driv_force); 
                   
        end
    end   
 
    %restarting the spacee loops to update the phi field
    for i=1:n
        for j=1:n
            phi(i,j) = phi(i,j) + update_phi(i,j);
        end
    end     
end    
    
phi_profile = phi(n/2,:);
%viewing the order parameter profile
plotres(phi_init,phi,phi_init_profile,phi_profile);

%imagesc(phi);
%colormap(gray)
%axis equal
%axis off




