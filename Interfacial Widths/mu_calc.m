%this function will take in the composition array and also the array index 
function[mu] = mu_calc(c,i,j)

%declaring the global variables
global dx
global dy
global H_b
global kappa
global sqdim


%calculating dfb_dc
dfb_dc = 2.0*H_b*c(i,j)*(1.0-c(i,j))*(1.0-2.0*c(i,j));

%calculating the laplacian
lapl_c = ((c(pbc(i+1,sqdim),j) - (2.0*c(i,j)) + c(pbc(i-1,sqdim),j))/(dx^2)) + ((c(i,pbc(j+1,sqdim)) - (2.0*c(i,j)) + c(i,pbc(j-1,sqdim)))/(dy^2));
mu = dfb_dc - kappa*lapl_c;

