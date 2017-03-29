function [w, output, Xs, sdev, Atrials] = sensitivity_AS(max_vals, min_vals, Nsamples, func, varargin)
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
    % 'w' contains the weight vector
    % 'output' contains the quantity of interest at each randomly selected
    %   parameter value
    % 'Xs' contains the normalized weights used to compute each randomly
    %   selected parameter
    % 'sdev' contains the standard deviations of the estimated means for
    %   each data point when using averaging
    % 'Atrials' contains the number of trials used for each data point when
    %   using averaging based on a tolerance
    
    % Ensure parameters are usuable
    parser = inputParser;
    addRequired(parser, 'max_vals', @isnumeric)
    addRequired(parser, 'min_vals', @isnumeric)
    addRequired(parser, 'Nsamples', @isnumeric)
    addRequired(parser, 'func')
    addParameter(parser, 'test_params', true(1, length(max_vals)), @islogical)
    addParameter(parser, 'Averaging', false, @islogical)
    addParameter(parser, 'Asamples', 10, @isnumeric)
    addParameter(parser, 'Atolerance', Inf, @isnumeric)
    
    parse(parser, max_vals, min_vals, Nsamples, func, varargin{:});
    test_params = parser.Results.test_params;
    Averaging = parser.Results.Averaging;
    Asamples = parser.Results.Asamples;
    Atolerance = parser.Results.Atolerance;
    
    output = zeros(Nsamples,1);   % Output of interest
    Nparams = sum(test_params);   % Number of parameters of interest
    sdev = zeros(Nsamples,1);     % Standard deviations from averaging
    Atrials = zeros(Nsamples,1);  % Number of trials used for averaging
    Xs = zeros(Nsamples,Nparams);       % To save the normalized paramters
    
    % Set base values for constant parameters
    params = max_vals;
    
    for jj = 1:Nsamples
        jj
        % Randomly select parameters from their spaces
        Xs(jj,:) = 2*rand(1, Nparams) - 1;
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
        sdev(jj) = std(out)/sqrt(ii);
        Atrials(jj) = ii-1;
        
    end
    % Solve the linear system
    A = [ones(Nsamples, 1), Xs];
    w = A\output;
    % Ignore the intercept
    w = w(2:end);
    % Normalize the weight vector
    w=(1/norm(w))*w;
end