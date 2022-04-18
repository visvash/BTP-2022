function[p]= plotres(phi_init,phi,phi_init_profile,phi_profile)


subplot(2,2,1)
imagesc(phi_init);colormap(gray);axis equal;axis off;title('initial','fontsize',18)

subplot(2,2,2)
imagesc(phi);colormap(gray);axis equal;axis off;title('final','fontsize',18)

subplot(2,2,3)
plot(phi_init_profile,'-o');xlabel('x','fontsize',18);ylabel('phi','fontsize',18)    

subplot(2,2,4)
plot(phi_profile,'-o');xlabel('x','fontsize',18);ylabel('phi','fontsize',18)    
