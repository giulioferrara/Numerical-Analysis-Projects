% ESERCIZIO 3, GRUPPO O 
% Camilla Mammoli, Giulio Ferrara, Giovanni Ufrasi

%% SCRIPT
% x-coordinates
x = linspace( -1, 4, 101 ); 

% y-coordinates
y = zeros( 101, 1 );
rng(1); eps = rand( 101, 1 );
y(:) = atan( x(:) ) + 0.15*( eps(:)-0.5 );

% Points used to build the polynomial
xk = x(1:2:end)';
yk = y(1:2:end);

% Defining the Vandermonde matrix
powers = 0:25;
A25 = xk.^powers;


% Loop on the three maximum degrees of the polynomial
for j = [25, 15, 6]

    % Submatrix of A25 nx(j+1)
    A = A25( :, 1:j+1 );

    % QR-factorization of A
    QR = householder(A);

    % Solving Ra = g to find polynomial coefficients
    R = triu(QR);
    g = find_g(QR, yk);

    % solving a = upper_triangular(R, g);
    a = g(1:j+1);

    for i = j+1:-1:1
        if R(i, i) == 0
            error('elemento diagonale nullo')
        end
        a(i) = a(i) - R(i, i+1 : j+1)* a(i+1 : j+1);
        a(i) = a(i)/ R(i,i);
    end


    % Flipping 'a' to have the coefficients in the correct order
    a = ( flip(a) )';

    % Plotting the correct function
    figure(j);
    
    xticks auto
    yticks auto
   
    v = linspace(-1, 4, 5000);
    w = atan(v);
    plot(v, w ,'k-', 'LineWidth', 1.7);
    hold on;
    grid on;
    
    % Plotting the approximated polynomial
    approx_pol = polyval(a, v);
    plot(v, approx_pol, 'b-', 'LineWidth', 4);
    hold on;
    
    % Even points
    xi = x(2:2:end)';
    yi = y(2:2:end);

    % Plotting odd and even points
    plot( xk, yk, "xmagenta", 'MarkerSize', 7, 'LineWidth', 1.7 );
    plot( xi, yi, "or", 'MarkerSize', 7, 'LineWidth', 1.7 );
    legend('atan(x)', 'approximate polynomial', 'odd points', ...
        'even points', 'Location', 'southeast', 'fontsize', 17);

   
    % Calculating RMS for odd and even points 
    oddval = polyval(a, xk);
    evenval = polyval(a, xi);

    RMSodd = sqrt( mean( (yk - oddval).^2) ) ;
    RMSeven = sqrt( mean( (yi - evenval).^2) );
  
    xlabel( ['RMS odd points: ', num2str(RMSodd), newline, newline, ...
         'RMS even points: ', num2str(RMSeven)], 'FontSize', 14.4 );
    
    
    
end

%% FUNCTION A
function [QR] = householder(A)
    
    % householder operates the QR factorization of A using 
    % the Householder method
    %
    % INPUTS
    % -A = mxn matrix (m>n)
    %
    % OUTPUTS
    % -QR = mxn matrix, containing the compressed form of Q and R

    QR = A;
    [m,n] = size(QR);

    for i = 1:n

        % Calculating alfa 
        alfa = norm( QR(i:m, i) );

        if alfa == 0
            error('A doesnt have maximum rank')
        end

        % Adjusting the sign of alfa
        if QR(i, i) >= 0
            alfa = -alfa;
        end
        
        % Vectorized update of QR
        v1 = QR(i, i) - alfa;
        QR(i, i) = alfa;

        QR(i+1:m, i) = QR(i+1:m, i) / v1;
        beta = -v1 / alfa;
        QR(i:m, i+1:n) = QR(i:m, i+1:n) - ( beta*[ 1; QR(i+1:m, i) ])*...
            ([1 QR(i+1:m, i)']*QR(i:m, i+1:n) );
    end    
end


%% FUNZIONE PUNTO B
function [g] = find_g(QR, b)
    
    % find_g calculates the array g = Q'b using a vectorized 
    % approach on QR, without explicitly calculating the matrix Q
    %
    % INPUTS
    % -QR = mxn matrix (m>n)
    % -b = nx1 array
    %
    % OUTPUTS
    % -g = nx1 array

    [m,n]=size(QR);

    % Initialization of g
    g = b;

    % Vectorized expression of g using the householder vectors in QR
    for i = 1:n

        % Householder vector
        v = [1; QR(i+1:m, i)];
        
        % Updating g
        beta = 2/ (v'*v);
        g(i:end) = g(i:end) - beta*v*( v'*g(i:end) );
    end
end

%%
%{
COMMENTI:
Si può osservare che all'aumentare del grado del
polinomio lo scarto quadratico medio tra quest'ultimo e i punti
pari aumenta, mentre per i punti dispari diminuisce. 
Questo è dovuto al fenomeno di overfitting: si perde l'andamento
generale dei punti e di questo risente l'RMS dei punti pari, i quali non
sono stati utilizzati per trovare il polinomio.
%}