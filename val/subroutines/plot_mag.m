function plot_mag(u, v)

uu = sqrt(u.^2+v.^2);
imagesc(uu);

axis equal off; 
drawnow;