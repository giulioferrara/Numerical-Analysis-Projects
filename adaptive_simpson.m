%ESERCIZIO 1, GRUPPO 'O' 
% Camilla Mammoli, Giulio Ferrara, Giovanni Ufrasi 

%% SCRIPT
% Integrand function
f = @(x) atan(x) - sqrt(x+1)/2;

% Interval
a = 0;
b = 8;




dfi=[3.09,23.8,104.5,366,1980,7385,27861,102700,242000,508900,1361000,3359000,5280000]...
    *pi*2./[0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20];
disp(dfi);



% Exact  value of the integral
format long;
true_val = integral(f, a, b, 'AbsTol', 1e-10, 'RelTol', 1e-10);
disp('Exact value of the integral');
disp(true_val);

% Parameters for the three figures
tols = [1e-5, 1e-6, 1e-7];
Ns = [64, 128, 256];

% Loop on the three configurations
for i = 1:3
    tol = tols(i);
    N = Ns(i);

    % Calculation of the integral with adaptive simpson
    [approx_val, nodes] = adaptive_simpson(f, a, b, tol, N);
    fprintf( 'Approximative value of the integral %d:\n', i);
    disp( approx_val );

    % Relative error
    rel_err = abs(approx_val - true_val)/abs(true_val);

    % Number of nodes used
    numnodes = length(nodes);

    % Plot
    figure(i);

    % Function
    x = linspace( a, b, 10000 ); 
    y = f(x);
    plot( x, y, 'b-','LineWidth', 5 );
    hold on;

    % Nodes
    xticks(nodes); 
    xticklabels ({});
    y_vals = f(nodes);
    plot( nodes, zeros(numnodes,1), 'ko', 'MarkerSize', 4 , 'MarkerFaceColor', 'y', 'LineWidth', 1.5 );
    plot( nodes, y_vals, 'bo', 'MarkerSize', 7, 'MarkerFaceColor', 'y', 'LineWidth', 2);
    
    grid on;

    % Title of each figure
    title(sprintf('tol = %.0e, N = %d', tol, N));

    % Text under the x axis
    xlabel(['Numero di nodi: ', num2str(numnodes),...
        '      Errore relativo: ',num2str(rel_err)]);
    
    % Text on the y axis
    ylabel('f(x)');

    % Legend
    legend('f(x)', 'nodes', 'f(nodes)', 'Location', 'southeast', 'FontSize', 17);

end


%% FUNZIONE PUNTO A 
function [I, nodes] = adaptive_simpson(f, a, b, tol, N)
    
    % adaptive_simpson approximates the integral of f between 
    % a and b. It operates by recursion using the auxiliary function
    % simpson_recursive
    %
    % INPUTS
    % -f = integrand function
    % -a = scalar, first end of integration interval
    % -b = scalar, last end of integration interval
    % -tol = scalar, tolerance for arrest criteria
    % -N = scalar, maximum number of intervals
    %
    % OUTPUTS
    % -I = scalar, approximate value of the integral
    % -nodes = array of used nodes
    
    
    h0 = b - a;
    fa = f(a);
    fb = f(b);
    c = (a + b) / 2;
    fc = f(c);
    
    format long;
    [I, nodes] = simpson_recursive(f, a, b, c, tol, N, h0, fa, fb, fc);

    % Adding the missing nodes
    nodes = [a, nodes, b]; 

    function [I, nodes] = simpson_recursive(f, a, b, c, tol, N, H, fa, fb, fc)
        
        % simpson_recursive uses the adaptive simpson method combined with
        % bisection to approximate the integral of f between a and b.
        %
        % INPUTS
        % -f, a, b, tol, N = same as in adaptive_simpson
        % -c = scalar, middle point between a and b
        % -H = b-a, scalar
        % -fa = f(a), scalar    fb = f(b), scalar    fc = f(c), scalar
        %
        % OUTPUTS
        % - I, nodes = same as in adaptive_simpson

        h = (b - a)/2;
        d = (a + c)/2;
        e = (c + b)/2;

        fd = f(d);
        fe = f(e);

        % Simpson on [a,b] and [a,c] + [c,b]
        I2 = (h/3) * (fa + 4*fc + fb);
        I4 = (h/6) * (fa + 4*fd + 2*fc + 4*fe + fb);
        

        if abs(I4 - I2) < 15 * tol      % Exit condtion
            I = I4;

            % Nodes used here
            nodes = [d c e];

        elseif (b - a) < H / N          % Error condition
            error(['[%f, %f] interval too short ' ...
                'compared to (b-a)/N'], a, b);

        else
            % Recursive call
           [Il, nodes_l] = simpson_recursive(f, a, c, d, ...
               tol/2, N, H, fa, fc, fd);

           [Ir, nodes_r] = simpson_recursive(f, c, b, e, ...
               tol/2, N, H, fc, fb, fe);

           % Total integral
           I = Il + Ir;

           % Sorted union of nodes
           nodes = [nodes_l, c, nodes_r];
        end
    end
end
%%
%{
COMMENTO:

Il criterio di arresto adottato nel metodo Simpson adattivo si dimostra efficace,
in quanto adatta dinamicamente la densità dei nodi in base alla difficoltà locale
della funzione integranda. Nei tratti più "curvi" o con più rapida variazione della funzione, il controllo
sull'errore (basato sulla differenza tra le due approssimazioni Simpson) forza
una suddivisione più fine, aumentando i nodi per mantenere la tolleranza richiesta.
Viceversa, dove la funzione è più regolare, il criterio riconosce che 
è possibile ottenere una buona approssimazione con pochi intervalli, evitando 
suddivisioni inutili. Questo garantisce un uso efficiente delle risorse computazionali,
mantenendo l’errore entro i limiti desiderati senza eccedere nel numero di nodi.
Questo comportamento evidenzia
l'efficacia del criterio di arresto adottato nella procedura di bisezione,
perché l'intervallo viene diviso solo quando è necessario.

Queste osservazioni valgono per tutte le combinazioni di N e tolleranza
considerate, e come previsto, all’aumentare di N o alla diminuzione della
tolleranza, i nodi si distribuiscono in maniera più fitta.
%}