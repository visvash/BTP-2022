% this is a code to simulate ferromagnetic domain growth
% we are going to decide the initial condititon
% all the simulation parameters are non-dimensional
function[]=Spinodal(n,tot_tsteps)

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

sqdim= n; %the grid dimensions are stored here  
dx=1; %the grid spacing along x
dy=1; %the grid spacing along y
dt=0.01; %the timestep size
H_b=1.0; %the parameter which sets the height of the activation barrier
kappa=1.0; %the gradient energy density coeffcient
M=1.0; %the mobility of the solutes 

%creating the initial condition 
c = 0.2 + 0.05*(2.0*randn(n) - 1);
for i=1:n
        for j=1:n
            if  (sqrt((i-(n/2))^2 + (j -(n/2))^2))<= 10
                c(i,j) = 1.0; 
            end
        end 
end

          
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
    e=(n/2);
    for i=(n/2):n
        if c(i,n/2)<0.5
            e=i;
            break
        end
    end
        m = c(e,n/2)-c(e-1,n/2);
        x = (0.5 - c(e,n/2)+ m*e)/m;
        tim(t)=t*0.01;
        radius(t) = (x-(n/2))^2;

end
%scatter(tim,radius)
m1=(c(e+1,n/2)-c(e-1,n/2))/2*dx;
m2=(c(e,n/2)-c(e-2,n/2))/2*dx;
x1= e - (c(e,n/2)/m1);
x2= (e-1) - (c(e-1,n/2)/m2);
y1=e+((1-c(e,n/2))/m1);
y2=(e-1)+((1-c(e-1,n/2))/m2);
xd=(x1+x2)/2;
yd=(y1+y2)/2;
Width = (xd-yd);
disp(H_b)
disp(kappa)
disp(Width)
%viewing the order parameter profile
%plotres(c_init, c);
%imagesc(phi);
%colormap(gray)
%axis equal
%axis off



