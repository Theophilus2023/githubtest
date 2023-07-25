function HopfBifurcation_sir()
    % Parameters
    beta = 0.25;     % Transmission rate
    gamma = 0.1;    % Recovery rate
    mu = 0.05;      % Natural birth/death rate
    nu = 0.005;     % Rate of loss of immunity
    R0 = beta/gamma;  % Basic reproduction number

    % Time span for simulation
    tspan = linspace(0, 200, 1000);

    % Initial conditions
    S0 = 0.9;   % Initial proportion of susceptible individuals
    I0 = 0.1;   % Initial proportion of infectious individuals
    R0 = 0.0;   % Initial proportion of recovered/immune individuals

    % ODE solver options
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

    % Solve the system of differential equations
    [t, sol] = ode45(@(t, y) sir_model(t, y, beta, gamma, mu, nu), tspan, [S0; I0; R0], options);

    % Extract S, I, and R components
    S = sol(:, 1);
    I = sol(:, 2);
    R = sol(:, 3);

    % Plot the results
    figure;
    plot(t, I, 'b-', 'LineWidth', 2);
    xlabel('Time (days)');
    ylabel('Proportion of Infectious Individuals');
    title('Hopf Bifurcation for Measles');
    grid on;
end

function dydt = sir_model(~, y, beta, gamma, mu, nu)
    S = y(1);
    I = y(2);
    R = y(3);

    % System of differential equations
    dSdt = mu - beta * S * I - mu * S;
    dIdt = beta * S * I - (gamma + mu + nu) * I;
    dRdt = gamma * I - mu * R + nu * I;

    dydt = [dSdt; dIdt; dRdt];
end
