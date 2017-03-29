function [evalues, U, output, outputplus, Xs, graddamp, sdev, Atrials] = Active_Subspaces(max_vals, min_vals, h, Nsamples, func, varargin)
    % INPUTS:
    % 'max_vals' is a vector containing the maximal values of the 
    %   parameters
    % 'min_vals' is a vector containing the minimal values of the
    %   parameters
    % 'h' is the finite differencing step size for the sensitivity metric
    % 'Nsamples' is the number of samples points for the Monte Carlo method
    % 'func' is the function of interestthat accepts a parameter vector and
    %   outputs a singular quantity dependent upon those parameters
    % 'test_params' is an optional vector containing the locations of the
    %   parameters for the sensitivity metric to be considered with respect
    %   to (non-constants)
    % 'Averaging' is a boolean that represents whether a number of samples
    %   should be taken for each data point (used for non-deterministic
    %   functions to reduce inaccuracy).  The default is to consider all
    %   parameters
    % 'Asamples' represents the (minimum) number of samples to take for 
    %   each data point.  The default is ten, and it will be considered 
    %   before the tolerance if a tolerance is chosen
    % 'Atolerance' represents the maximum standard deviation allowed for
    %   the expected value of a data point
    %
    % OUTPUTS:
    % 'evalues' contains the eigenvalues of the computational gradient
    %   matrix in descending order.  Eigenvalues are computed using the SVD
    %   command
    % 'U' contains the associated eigenvectors - U(:, 1) contains the
    %   eigenvector corresponding to the largest eigenvalue
    % 'output' contains the quantity of interest at each randomly selected
    %   parameter value
    % 'outputplus' contains the quantity of interest at each perturbation
    %   of the randomly selected parameter values
    % 'Xs' contains the normalized weights used to compute each randomly
    %   selected parameter
    % 'graddamp' is the computational gradient matrix
    % 'sdev' contains the standard deviations of the estimated means for
    %   each data point when using averaging
    % 'Atrials' contains the number of trials used for each data point when
    %   using averaging based on a tolerance
    
    % Ensure parameters are usuable
    parser = inputParser;
    addRequired(parser, 'max_vals', @isnumeric)
    addRequired(parser, 'min_vals', @isnumeric)
    addRequired(parser, 'h', @isnumeric)
    addRequired(parser, 'Nsamples', @isnumeric)
    addRequired(parser, 'func')
    addParameter(parser, 'test_params', true(1, length(max_vals)), @islogical)
    addParameter(parser, 'Averaging', false, @islogical)
    addParameter(parser, 'Asamples', 10, @isnumeric)
    addParameter(parser, 'Atolerance', Inf, @isnumeric)
    
    parse(parser, max_vals, min_vals, h, Nsamples, func, varargin{:});
    test_params = parser.Results.test_params;
    Averaging = parser.Results.Averaging;
    Asamples = parser.Results.Asamples;
    Atolerance = parser.Results.Atolerance;
    
    output = zeros(Nsamples,1);            % Output of interest
    Nparams = sum(test_params);         % Number of parameters of interest
    outputplus = zeros(Nparams,Nsamples);  % Perturbed output of interest
    sdev = zeros(Nparams+1, Nsamples);     % Standard deviations from averaging
    Atrials = zeros(Nparams+1, Nsamples);  % Number of trials used for averaging
    graddamp = zeros(Nparams,Nsamples);    % Gradient of output of interest 
    Xs = zeros(Nsamples,Nparams);          % To save the normalized paramters
    evalues = zeros(Nparams,1);            % Eigenvalues of the C matrix
    I = eye(Nparams);                      % Identity matrix
    
    % Set base values for constant parameters
    params = max_vals;
    
    for jj = 1:Nsamples
        jj
        % Randomly select parameters from their spaces
        Xs(jj, :) = 2*rand(1, Nparams) - 1;
        params(test_params) = 1/2*(diag(max_vals(test_params) - min_vals(test_params))*Xs(jj, :)' + (max_vals(test_params) + min_vals(test_params))');
        
        % Average the quantities of interest computed for each parameter
        %   value
        not_done = true;
        ii = 1;
        while not_done
            out(ii) = func(params);
            not_done = Averaging&((ii<=Asamples) || (std(out)/sqrt(ii) > Atolerance));
            ii = ii + 1;
        end
        
        output(jj) = mean(out);
        sdev(1, jj) = std(out)/sqrt(ii);
        Atrials(1, jj) = ii-1;
        
        for kk = 1:Nparams
            % Perturb the previously selected parameters
            xplus = Xs(jj, :)' + h*I(:, kk);
            paramsplus(test_params) = 1/2*(diag(max_vals(test_params) - min_vals(test_params))*xplus + (max_vals(test_params) + min_vals(test_params))');

            not_done = true;
            ii = 1;
            while not_done
                outplus(ii) = func(paramsplus);
                not_done = Averaging&((ii<=Asamples) || (std(outplus)/sqrt(ii) > Atolerance));
                ii = ii+1;
            end
            
            outputplus(kk, jj) = mean(outplus);
            sdev(kk+1, jj) = std(outplus)/sqrt(ii);
            Atrials(kk+1, jj) = ii-1;
        end
        % Store the results in the gradient matrix
        graddamp(:, jj) = (outputplus(:,jj) - output(jj))/h;
    end
    % Correct for the Monte Carlo estimation
    graddamp = 1/sqrt(Nsamples)*graddamp;
    % Compute the SVD decomposition of the computational gradient matrix
    [U, S, V] = svd(graddamp);
    evalues = diag(S.^2);
end