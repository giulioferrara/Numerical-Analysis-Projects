% ESERCIZIO 4, GRUPPO O 
% Camilla Mammoli, Giulio Ferrara, Giovanni Ufrasi

%% SCRIPT
% Creating sparse matrix A0
rng(1);
n = 1000;
A0 = speye(n) + sprand(n, n, 1/n);

% Defining A
A = A0'*A0;

% Exact solution
c = 1:1000;

% Known term
b = A*c';

% Toletance for arrest criteria
tol = 10^(-6);

% Conditioning number
cond_A = condest(A,2);

% Printing results
fid = fopen('risultati.txt', 'w');
fprintf(fid, 'Numero di condizionamento: %.2e\n', cond_A);
fprintf(fid, '<><><><><><><><><><><><><><><><><><>\n');

% Loop for the three values of omega
for w=[0.75, 1, 1.35]

    % Approximated solution and number of iterations
    [x,count] = relSOR(A, b, w, tol, n);

    % Relative error
    err_relativo = norm( A*x-b )/ norm(b);
    fid = fopen('risultati.txt', 'a');
    fprintf(fid, 'omega = %.2f\n', w);
    fprintf(fid, 'numero di iterazioni = %d\n', count);
    fprintf(fid, 'Errore relativo: %.2e\n', err_relativo);
    fprintf(fid, '-------------------------------\n');
    fclose(fid);
end

%% FUNCTION A
function [xk, count] = relSOR(A, b, w, tol, maxiter)
    
    % relSOR solves the linear system Ax=b with the iterative method
    % of Gauss-Siedel with relaxation parameter w
    %
    % INPUTS
    % -A = nxn matrix
    % -b = nx1 column array, known term
    % -w = scalar, in (0, 2)
    % -tol = scalar, tolerance for arrest criteria
    % -maxiter = scalar integer, maximum number of iterations
    %
    % OUTPUTS
    % -xk = nx1 array, approximated solution
    % -count = number of executed iterations

    n = size(A,1);

    % defining D, L, M
    D = spdiags(diag(A), 0, n, n);
    L = tril(A, -1);
    M = D + w*L;

    % Changing A and b to solve wA*x=wb
    b = w*b;
    A = w*A;

    % Starting point
    xk = zeros(n, 1);

    % First residual
    r = b;

    for count = 1:maxiter

        % Finding difference between consecutive values of xk
        u = r;
        for i = 1:n
            u(i) = u(i)- M(i, 1:i-1)*u(1:i-1);
            u(i) = u(i)/M(i, i);
        end

        % Updating xk
        xk = xk + u;

        % Current residual and arrest criteria 
        r = b - A*xk;

        if norm(r)/ norm(b)<=tol
            break
        end
        
    end
    
    if count == maxiter
        error('max number of iterations reached, convergence not guaranteed')
    end

end


%%
%{ 
COMMENTI:
Con w = 0.75, il metodo converge più lentamente: sottorilassare 
rallenta la convergenza.
Con w = 1, si ottiene il metodo di Gauss-Seidel classico e già si 
osserva un notevole miglioramento: la convergenza è più veloce 
rispetto a ω = 0.75.
Con w = 1.35, si ha un ulteriore miglioramento in efficienza: il numero 
di iterazioni diminuisce drasticamente, mostrando come un buon
sovrarilassamento può accelerare fortemente la convergenza.
L'errore relativo è dello stesso ordine per tutti e tre gli omega, cambiano
solo i decimali della mantissa. Nessuno dei tre w compromette la
convergenza. 
Con delle prove abbiamo scoperto che molto probabilmente l'w migliore per
la nostra matrice è circa 1.79
La convergenza in realtà è assicurata anche dal fatto che A è sdp (e
0<w<2)infatti:
x^T * A * x = x^T * A0^T * A0 * x = (A0*x)^T (A0*x) = norm(A0*x)^2 > 0

perché  A0*x è diverso da 0  se lo è x. Infatti abbiamo verificato che A0 è
non singolare controllando i suoi autovalori tramite il comando eig(A0).

La velocità di convergenza migliora anche dal fatto che A si avviccina ad essere a
diagonale dominante e il metodo SOR funziona meglio per queste matrici.
%}

