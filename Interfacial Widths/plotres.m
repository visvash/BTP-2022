function[p]= plotres(c_init,c)


subplot(1,2,1)
imagesc(c_init);colormap(gray);axis equal;axis off;title('initial','fontsize',18)

subplot(1,2,2)
imagesc(c);colormap(gray);axis equal;axis off;title('final','fontsize',18)

%subplot(2,2,3)
%plot(phi_init_profile,'-o');xlabel('x','fontsize',18);ylabel('phi','fontsize',18)    

%subplot(2,2,4)
%plot(phi_profile,'-o');xlabel('x','fontsize',18);ylabel('phi','fontsize',18)    
