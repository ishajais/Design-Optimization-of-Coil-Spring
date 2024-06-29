function minimumMassSpringDesign()
    Q = 5;
    P = 1000;
    rho = 7850;
    G = 2e4; 
    delta = 0.5; 
    wo = 10; 
    tau_o = 250e6; 
    dmin = 0.1; 
    dmax = 1; 
    Dmin = 0.5; 
    Dmax = 2; 
    Nmin = 5; 
    Nmax = 10; 
    lb = [Nmin, Dmin, dmin];
    ub = [Nmax, Dmax, dmax];
    x0 = [Nmin, Dmin, dmin];
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', ...
                           'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true);

    [x, fval] = fmincon(@(x) objectiveFunction(x, Q, rho), x0, [], [], [], [], lb, ub, ...
                        @(x) constraints(x, P, G, tau_o, delta, wo, rho), options);
    fprintf('Optimal values: N = %.2f, D = %.2f, d = %.2f\n', x(1), x(2), x(3));
    fprintf('Minimum mass: %.4f kg\n', fval);
end

function [mass, grad] = objectiveFunction(x, Q, rho)
    N = x(1);
    D = x(2);
    d = x(3);
    mass = 0.25 * (N + Q) * pi^2 * D * d^2 * rho;
    grad = [
        0.25 * pi^2 * D * d^2 * rho;
        0.25 * pi^2 * (N + Q) * d^2 * rho;
        0.5 * pi^2 * (N + Q) * D * d * rho
    ];
end

function [c, ceq, c_grad, ceq_grad] = constraints(x, P, G, tau_o, delta, wo, rho)
    N = x(1);
    D = x(2);
    d = x(3);
    c = [
        -(8 * P * D^3 * N) / (d^4 * G) + delta; 
        (8 * P * D) / (pi * d^3) * ((4 * D - d) / (4 * (D - d)) + 0.615 * d / D) - tau_o; 
        -(d / (2 * pi * N * D^2) * sqrt(G / (2 * rho))) + wo; 
    ];
    c_grad = [
        -8 * P * D^3 / (d^4 * G), -24 * P * D^2 * N / (d^4 * G), 32 * P * D^3 * N / (d^5 * G);
        0, 8 * P / (pi * d^3) * (4 - 1 / (1 - d / (4 * D)) - 0.615 * d / D^2), ...
        -24 * P * D / (pi * d^4) * ((4 * D - d) / (4 * (D - d)) + 0.615 * d / D);
        0, 2 * d * sqrt(G / (2 * rho)) / (pi * N * D^3), -sqrt(G / (2 * rho)) / (2 * pi * N * D^2) ...
    ];
    ceq = [];
    ceq_grad = [];
end
