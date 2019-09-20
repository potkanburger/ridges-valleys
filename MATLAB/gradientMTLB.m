ITEMP = imread('smooth_lena.jpg');
I = rgb2gray(ITEMP);

[Gx,Gy] = imgradientxy(I, 'central');
figure(1)
imshowpair(Gx,Gy,'montage')
title('Directional Gradients Gx and Gy, Using Central Difference')
saveas(gcf,'lena_gradMATLAB.png')

pyXgrad = imread('x_gradient_smooth_lena.jpg');
pyYgrad = imread('y_gradient_smooth_lena.jpg');
figure(2)

imshowpair(pyXgrad,pyYgrad,'montage')
title('Directional Gradients Gx and Gy, Using Central Difference, From Python')
saveas(gcf,'lena_gradPYT.png')