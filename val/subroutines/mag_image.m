function ii = mag_image(u, v)

uu = sqrt(u.^2+v.^2);
ii = image(uu);