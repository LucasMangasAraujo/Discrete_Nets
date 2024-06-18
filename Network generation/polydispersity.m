%% polydispersity
% -------------------------------------------------------------------------
% This function creates the polydisperse twin of the generated network
% according to the input distribution. The function also generates a figure
% comparing the empirical distribution with the analytical one.
% 
% OBS: Current version of the function only supports the log-normal distribution. 
% However, other distributions can be readly added by simply addig a elif 
% statement.
% 
% distribution: String identifying the distribution that will be used
% par: Cell with the relavant parameters of the distribution
% Nbonds: number of chains in the network
% 
% The function returns
% 1) Array containing the chain lengths: chain_lengths
% -------------------------------------------------------------------------
function chain_lengths = polydispersity(distribution, par, Nbonds)
    if distribution == "Log"
        mu_N = par{1}; s_N = par{2}; % Mean and deviation of the chain length distribution

        % Calculate the parameters of the log-normal distribution
        mu = log((mu_N^2)/sqrt(mu_N^2 +  s_N^2));
        s = sqrt(log(1 + (s_N^2/mu_N^2)^2));

        % Draw Nbonds chain lentghs from the distribution
        chain_lengths = lognrnd(mu, s, [Nbonds, 1]);
        x = linspace(0.9 * min(chain_lengths),  1.1 * max(chain_lengths), 1000);
        analytic = lognpdf(x, mu, s);

    end

    % Compare analytic and empirical distributions
    figure(2)
    h = histogram(chain_lengths, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold on
    p = plot(x, analytic, 'b-', 'LineWidth', 2);
    xlabel('$N$', 'Interpreter', 'latex', FontSize=35);
    ylabel('$p_N$', 'Interpreter', 'latex', FontSize=35);
    legend([h p], {'Empirical', 'Analytic'}, 'Location', 'best', Interpreter='latex');
    set(gca, 'FontSize', 25);
    set(gcf,'Color','w');
    axis square
    hold off
end