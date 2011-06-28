function [] = EM_step_plot( mu, SIGMA, step, corners )

if any([1 2 5 10 17 26 34 37 38 39 40 41] == step) % only plot these steps

    res = 100;
    
    disp(['plotting step ', num2str(step)])

    plot_gauss_mix(mu, SIGMA, corners, 'res', res, 'kernel_indices', [16 334])
    
    title(['two selecter kernels after EM step ', num2str(step)])
    
    hold on
    
%    scatter(mu(:,1),mu(:,2),'.w')
    
    hold off
    

    axis off
    axis equal
    set(gcf, 'PaperPositionMode', 'auto', ...
             'units','normalized', ...
             'outerposition',[0 0 1 1]) % maximize before saving

end

end

