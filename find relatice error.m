

%c
for n = [10, 100, 1000]
    % initialize 
    b = ones(n,1);
    A = ones(n);
    soln_exact = zeros(n,1);
    for ii = 1 : n
        A(ii,ii)=n;
        soln_exact(ii) = 1/(2*n-1);
    end
    % get the returned data from Jacobi and GSeidel
    [xj,itj] = Jacobi(A,b);
    [xg,itg] = GSeidel(A,b);
    % find the relative error of the 
    errorj = (norm (xj-soln_exact)/norm(soln_exact));
    errorg = (norm (xg-soln_exact)/norm(soln_exact));
    fprintf(['n: %d, errorj: %d, itj: %d, errorg: %d, itg: %d,\n'],n,errorj,itj,errorg,itg);
end

%a

% Jacobi function
% Input :a square n*n matrix A and an n*1 vector b. 
% Return the iterations and a value x.
function [x, iteration] = Jacobi (A,b)
    % initialize
    n = size(b); % get n 
    x_old = zeros(n); % initial guess a zero vector
    error = 10^(-6); % tolerance
    Max_iter = 1000; % Maximum iterations
    iteration = 1; % initialize iteration

    while iteration < Max_iter
        % initialize x_new
        x_new = zeros(n);
        
        % get x_new
        for ii = 1: n

            %get the sum of aij*xj
            sum = 0;
            for jj = 1 : n
                if jj ~= ii
                    sum = sum +A(ii,jj) * x_old(jj);
                end
            end

            x_new(ii)=(1/A(ii,ii))*(b(ii)-sum);
        end
        
        %check tolerance
        if norm(x_new - x_old) < error
            break
        end

        iteration = iteration + 1;
        x_old = x_new;
    end
    % return the x value
    x = x_new;
end


%b

% GSeidel function
% Input :a square n*n matrix A and an n*1 vector b. 
% Return the iterations and a value x.
function [x, iteration] = GSeidel(A,b)
    % initialize
    n = size(b); % get n 
    x_old = zeros(n); % initial guess a zero vector
    error = 10^(-6); % tolerance
    Max_iter = 1000; % Maximum iterations
    iteration = 1; % initialize iteration
    x_new = x_old; % initial x_new

    while iteration < Max_iter
        % get x_new
        for ii = 1: n

            %get the sum of aij*xj new
            sum_new = 0;
            for jj = 1 : ii-1
                sum_new = sum_new +A(ii,jj) * x_new(jj);
            end
            
            %get the sum of aij*xj old
            sum_old = 0;
            for j = ii+1 : n
                sum_old = sum_old + A(ii,j) * x_old(j);
            end
             
            x_new(ii)=(1/A(ii,ii))*(b(ii)-sum_new - sum_old);
        end
        
        %check tolerance
        if norm(x_new - x_old) < error
            break
        end

        iteration = iteration + 1;
        x_old = x_new;
    end
    % return the x value
    x = x_new;
end
      
