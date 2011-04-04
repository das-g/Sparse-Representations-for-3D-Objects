function [ ] = kernel_step_plot( mu, SIGMA, step, corners)
%KERNEL_STEP_PLOT Summary of this function goes here
%   Detailed explanation goes here

plot_steps = [1 10 40]; % Steps we want to plot
observe_idx = 45; % index of the kernel to be observed

if any(step == plot_steps) % plot the current step?
    plot_gauss_mix(mu, SIGMA, corners, ...
                   'res', 400, ...
                   'kernel_indices', observe_idx)
    hold on
    scatter(mu(:, 1), mu(:, 2))
    scatter(mu(observe_idx, 1), mu(observe_idx, 2), 'red')
    
    saveas(gcf, ... current figure
           ['kernel_EM-step' int2str(step) '.fig'])
    close(gcf)
end

end

