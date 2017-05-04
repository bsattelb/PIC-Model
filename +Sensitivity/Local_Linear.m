function [evalues, U, output, Xs, Xs2, graddamp, sdev, Atrials] = Local_Linear_Model(max_vals, min_vals, Nsamples, Nsamples2, p, func, varargin)
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
    
    output = zeros(Nsamples,1);     % Output of interest
    Nparams = sum(test_params);     % Number of parameters of interest
    sdev = zeros(Nsamples,1);       % Standard deviations from averaging
    Atrials = zeros(Nsamples,1);    % Number of trials used for averaging
    Xs = zeros(Nsamples,Nparams);   % To save the normalized parameters used in calculating the QOI
    Xs2 = zeros(Nsamples2,Nparams); % To save the normalized parameters used for local linear
    b = zeros(Nsamples2,Nparams);   % non-constant coefficients of the local linear regression models
    
    % Set base values for constant parameters
    params = max_vals;
	
	% Convert from normalized parameters to values used in the function of interest
	convertParams = @(x) 1/2*(diag(max_vals(test_params) - min_vals(test_params))*x + (max_vals(test_params) + min_vals(test_params))');
    
    for jj = 1:Nsamples
        jj
        % Randomly select parameters from their spaces
        Xs(jj,:) = 2*rand(1, Nparams) - 1;
        params(test_params) = convertParams(Xs(jj, :)');
        
        % Compute mean value of quantities of interest for each
		%   given parameter value
        not_done = true;
        while not_done
			Atrials(jj) = Atrials(jj) + 1;
            out(Atrials(jj)) = func(params);
            not_done = Averaging&((Atrials(jj)<=Asamples) || (std(out)/sqrt(Atrials(jj)) > Atolerance));
        end
        
        output(jj) = mean(out);
        sdev(jj) = std(out)/sqrt(Atrials(jj));
        
    end
    
    for jj = 1:Nsamples2
        % Randomly select parameters
        Xs2(jj,:) = 2*rand(1, Nparams) - 1;
        
        % Find the p nearest elements of Xs
        dif = bsxfun(@minus, Xs2(jj,:), Xs);
        dist = sqrt(sum(dif.^2, 2));
        [~, Ix] = sort(dist, 'ascend');
        NstNghbrXS = Xs(Ix(1:p),:);
        NstNghbrQOI = output(Ix(1:p));
        
        % Find the local linear regression model
        bBig = [ones(p,1), NstNghbrXS]\NstNghbrQOI;
        
        % Truncate to ignore the constant coefficient
        b(jj,:) = bBig(2:end);
    end
    
    graddamp = 1/sqrt(Nsamples2)*(b'*b);
    % Compute the SVD decomposition of the computational gradient matrix
    [U, S, V] = svd(graddamp);
    evalues = diag(S.^2);
end