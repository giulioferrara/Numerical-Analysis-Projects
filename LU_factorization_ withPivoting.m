%ESERCIZIO 2, GRUPPO 'O' 
% Camilla Mammoli, Giulio Ferrara, Giovanni Ufrasi 

%% SCRIPT
fid = fopen('report.txt', 'w');

for n=10:2:20
    fid = fopen('report.txt', 'a');

    % Assembling of c and r, column and row vectors
    c = nthroot(3, n:-1:1);
    r = 3.^(1:n);

    % Exact solution
    xhat=(n:-1:1)';
    
    % Creation of the Hankel matrix of dimension nxn
    H = hankel(c, r);

    % Known term of the linear sistem
    b = H*xhat;

    % Approx solution with partial pivoting
    x0 = lin_solve(H,b,0);

    % Approx solution with total pivoing
    x1 = lin_solve(H,b,1);

    % Conditioning number of the Hankel matrix
    cond_H = cond(H,2);

    % Relative errors for partial and total pivoting
    partialp_err = norm(xhat-x0,2) / norm(xhat,2);
    totalp_err = norm(xhat-x1,2) / norm(xhat,2);

    % Printing results on the report file
    fprintf (fid, 'n = %d\n', n);
    fprintf (fid, 'Conditioning number: %.2e\n', cond_H);
    fprintf (fid, 'Relative error for partial pivoting: %.2e\n', ...
        partialp_err);
    fprintf (fid, 'Relative error for total pivoting: %.2e\n', ...
        totalp_err);
    fprintf (fid, '-------------------------------\n');

    fclose (fid);
    
end



%% FUNCTION A
function [LU, p, q] = total_pivoting(A)

    % total_pivoting operates the LU factorization on the matrix A 
    % with total pivoting
    %
    % INPUTS
    % -A = nxn matrix
    %
    % OUTPUTS
    % -LU = nxn matrix containing the compressed form of L and U,
    % such as L*U=A
    % -p = 1xn row vector, permutation of the rows of the matrix A
    % -q = 1xn row vector, permutation of the columns of the matrix A

    n = size(A,2);
    assert( size(A,1) == n, 'non-square matrix!')
    
    LU = A;
    
    % Initialisation of the permutation vectors p and q
    p = 1:n;
    q = 1:n;

    for k = 1:n-1
        
        % Finding the absolute maximum of the kxk submatrix 'sub'
        % and its indexes i,j in the matrix LU
        sub = LU(k:n,k:n);
     
        s = size(sub,1);

        %finding maximum of the flattened by columns matrix 
        [mu,ind] = max(abs (sub (:) ) );

        % discovering the indeces
        [i,j]=ind2sub([s,s],ind);
        
        i = i+k-1;
        j = j+k-1;

        if mu == 0
            error('The matrix is singular')
        end


        % Permutation of rows and columns of LU if needed 
        if i > k
            p( [i k] ) = p( [k i] );
            LU( [i k], :) = LU( [k i], :);
        end

        if j > k
            q( [j k] ) = q( [k j] );
            LU( :, [j k] ) = LU( :, [k j] );
        end
        
        % Vectorized update of LU 
        LU( k+1:n, k ) = LU( k+1:n, k )./ LU( k, k );
        LU( k+1:n, k+1:n ) = LU( k+1:n, k+1:n )- ... 
        LU( k+1:n, k )*LU( k, k+1:n );
    end

    if LU(n,n)==0
        warning ('matrice fattorizzabile LU ma singolare')
    end


end




%% FUNCTION B
function [x, LU, p, q]= lin_solve(A, b, F)

   % lin_solve solves the linear system Ax=b using LU factorization of A
   % with partial (F=0) or total (F=1) pivoting
   %
   % INPUTS
   % -A = nxn matrix
   % -b = nx1 column vector, known term of the system
   % -F = scalar (0 or 1)
   %
   % OUTPUTS
   % -x = nx1 column vector, solution of the system
   % -LU = nxn matrix, compressed form of the factorization
   % -p = 1xn row vector, permutation vector
   % -q = 1xn row vector, permutation vector

   n = size(A,2);

    switch F
        case 0 %partial pivoting

            assert( size(A,1) == n, 'matrice non quadrata!')
            LU = A;

            % Initialization of permutation vector of the rows of LU
            p = 1:n;

            for k = 1:n-1

                % Finding absolute maximum value of the k-th column
                % of the kxk submatrix LU(k:n, k:n) and its index ind
                [mu, ind] = max( abs( LU( k:n, k ) ) );
                ind=ind+k-1;

                if mu == 0
                    error('matrice singolare')
                end
                
                % Permuting the rows of LU if needed
                if ind > k
                    p( [ind k] ) = p( [k ind] );
                    LU( [ind k], :) = LU( [k ind], :);
                end

                % Vectorized update of LU
                LU( k+1:n, k ) = LU( k+1:n, k )/ LU( k, k );
                LU( k+1:n, k+1:n ) = LU( k+1:n, k+1:n )-...
                LU( k+1:n, k )*LU( k, k+1:n );
            end
            q=1:n;
            
        case 1 %total pivoting

            [LU,p,q] = total_pivoting(A);
    end 
    
    % Finding L and U
    L = eye( size(LU, 1) )+tril(LU, -1);
    U = triu( LU );

    % Solving Ly=b(p)
    y = b(p);

    for i = 1:n

        if L(i, i) == 0
            error('null diagonal entry')
        end

        y(i) = y(i)- L(i, 1:i-1)*y(1:i-1);
        y(i) = y(i)/ L(i, i);
    end


    % Solving Ux=y
    x = y;
    for i = n:-1:1

        if U(i,i)==0
            error('null diagonal entry')
        end

        x(i) = x(i)- U(i, i+1:n)*x(i+1:n);
        x(i) = x(i)/ U(i, i);
    end


    % Permuting the solution x in case of total pivoting
    if F==1
        % qt is the inverse permutation of q
        qt(q)=1:n;
        x=x(qt);
        
    end

end

%% 
%{
COMMENTI:

Dal report risulta che il pivoting totale garantisce una precisione superiore
rispetto a quello parziale. Tuttavia, questo miglioramento non è particolarmente
rilevante rispetto all’ordine di grandezza complessivo degli errori. 
Di conseguenza, se si è disposti a tollerare una leggera perdita in accuratezza,
il pivoting parziale è da preferire per via del costo computazionale inferiore.
Al termine dell’algoritmo con pivoting totale, è importante ricordare che,
a differenza del pivoting parziale, durante le operazioni possono essere 
scambiate anche le colonne della matrice. Per questo motivo, 
la soluzione ottenuta x deve essere riorganizzata tramite una permutazione
inversa qt del vettore q che tiene traccia degli scambi tra colonne: 
quindi si deve porre x(qt) affinché le componenti della soluzione 
siano riportate nell’ordine corretto.
%}