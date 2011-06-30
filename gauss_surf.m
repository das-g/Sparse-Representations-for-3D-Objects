% Just some illustrative plots for the beamer presentation.

a = linspace(-1,1,30);
[X Y] = meshgrid(a,a);

figure
surf(X,Y,X)
set(gca, 'XTick', [0], ...
         'YTick', [0], ...
         'ZTick', [0])
saveas(gcf, 'surf_plot_linear', 'epsc')

figure
Z = gaussmf(sqrt(X.^2+Y.^2), [0.3 0]);
surf(X,Y,Z)
axis off
saveas(gcf, 'surf_plot_gauss', 'epsc')