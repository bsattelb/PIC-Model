function movie_run( DT, NT, NG, N, distribution, params, saveFrameNum, movieFrameNum)


    % Ensure that computational parameters are set up correctly
    parser = inputParser;
    this_dir = strrep(which(mfilename('fullpath')), [mfilename '.m'],'');
    distribution_files = dir([this_dir 'Initilization\*.m']);
    possible_distributions = regexprep(regexprep({distribution_files.name}, '.[^.]*$', ''), '^[^_]*_', '');
    addRequired(parser, 'DT', @isnumeric)
    addRequired(parser, 'NT', @isnumeric)
    addRequired(parser, 'NG', @isnumeric)
    addRequired(parser, 'N', @isnumeric)
    addRequired(parser, 'distribution', @(x) any(validatestring(x, possible_distributions)));

    parse(parser, DT, NT, NG, N, distribution);

    L = params(1);
    NGp = max(64, NG/2^3);
    dx = L/NG;

    curDir = pwd;
    cd([this_dir 'Initilization']);
    temp = ['initilization_' distribution];
    [xp, vp, rho_back, Q, QM] = feval(temp, params, N);
    cd(curDir);

    p=1:N;
    p=[p p];
    un=ones(NG-1,1);
    Poisson=spdiags([un -2*un un],[-1 0 1],NG-1,NG-1);
    time = 0:DT:NT*DT;

    histTotEnergy = zeros(1, NT);
    histKinEnergy = zeros(1, NT);
    histPotEnergy = zeros(1, NT);
    histMomentum = zeros(1, NT);
    histSpectrum = [];
    histEfield = [];
    histEL2 = zeros(1, NT);
    
    if ~exist(fullfile(cd, 'images'), 'dir')
        mkdir('images');
    end
    
    if ~exist('movies', 'dir')
        mkdir('movies');
    end
    
    f = figure(1);
    set(gcf, 'Visible', 'off');
    Title = {['Landau Damping']; ['dt = ' num2str(DT) ', N = ' num2str(N, '%g') ', dx = ' num2str(dx)]; ['T = ' num2str(DT*0)]};
    phaseSpace(xp, vp, NGp, Title);
    
    saveas(f, 'images/initialPhaseSpace.png');
    
    video = VideoWriter(['movies/phaseSpaceBehavior_' distribution '.mp4'], 'MPEG-4');
    video.FrameRate = 4;
    open(video);
    writeVideo(video, getframe(f));
    pause(0.01);
    
    
    for it=1:NT

        % update particle position xp
        xp=xp+vp*DT;  
        % apply the boundary conditions
        out=(xp<0); xp(out)=xp(out)+L;
        out=(xp>=L);xp(out)=xp(out)-L;

        % interpolation from particle to grid
        g1=floor(xp/dx-.5)+1;
        g=[g1;g1+1];
        fraz1=1-abs(xp/dx-g1+.5);
        fraz=[fraz1;1-fraz1]	;
        out=(g<1);g(out)=g(out)+NG;
        out=(g>NG);g(out)=g(out)-NG;
        mat=sparse(p,g,fraz,N,NG);
        % calculate the charge density
        rho=full((Q/dx)*sum(mat))'+rho_back;

        % calculate the electrostatic potential'
        Phi=[Poisson\(-rho(1:NG-1)*dx^2); 0];

        % calculate the electric field from the electrostatic potential
        Eg=([Phi(NG); Phi(1:NG-1)]-[Phi(2:NG);Phi(1)])/(2*dx);

        % interpolation grid -> particle and velocity update
        vp=vp+mat*QM*Eg*DT;

        % Calculate different kind of energies
        Ekin   = 0.5*abs(Q)*sum(vp.^2);
        % Potential energy
        Efield = 0.5*sum(Eg.^2)*dx;
        % Total energy
        Etot =  Ekin + Efield ;
        %L2 norm
        EL2 = sqrt(2*Efield);

        % Total momentum
        histMomentum(it) = sum(abs(Q)*vp);

        % history of various energies
        histKinEnergy(it) = Ekin;
        histPotEnergy(it) = Efield;
        histTotEnergy(it) = Etot;
        histEfield = [histEfield Eg];
        histEL2(it) = EL2;

        % take the Fourier transform of the electric field on the grid
        NFFT = 2^nextpow2(length(Eg)); % Next power of 2 from length of Eg
        Y = fft(Eg,NFFT)/length(Eg);
        histSpectrum = [histSpectrum 2*abs(Y(1:NFFT/2+1))];

        if mod(movieFrameNum*it*DT, 1) == 0
            Title = {['Landau Damping']; ['dt = ' num2str(DT) ', N = ' num2str(N, '%g') ', dx = ' num2str(dx)]; ['T = ' num2str(DT*it)]};
            phaseSpace(xp, vp, NGp, Title);
            writeVideo(video, getframe(f));
            pause(0.01);
            if mod(saveFrameNum*it*DT, 1) == 0
                saveas(f, ['images/T' num2str(it*DT), 'PhaseSpace.png']);
            end
        end
    end

    
    Title = {['Landau Damping']; ['dt = ' num2str(DT) ', N = ' num2str(N, '%g') ', dx = ' num2str(dx)]; ['T = ' num2str(NT*DT)]};
    phaseSpace(xp, vp, NGp, Title);
    saveas(gcf, 'images/finalPhaseSpace.png');
    
    close(video);

%     curDir = pwd;
%     cd([this_dir 'QOI_calc']);
%     temp = ['calc_' distribution];
%     damprate = feval(temp, time, histEnergy);
%     cd(curDir);
    
    semilogy(time(2:end), histSpectrum(2, :));
    title({[distribution ' Second Fourier Mode']; ['dt = ' num2str(DT) ', N = ' num2str(N, '%g') ', dx = ' num2str(dx)]});
    saveas(f, 'images/fourierDamping.png')
    
    semilogy(time(2:end), histEL2);
    title({[distribution ' ||E||_{L2}']; ['dt = ' num2str(DT) ', N = ' num2str(N, '%g') ', dx = ' num2str(dx)]});
    saveas(f, 'images/fieldDamping.png')
    
    semilogy(time(2:end), histKinEnergy);
    title({[distribution ' Kinetic Energy']; ['dt = ' num2str(DT) ', N = ' num2str(N, '%g') ', dx = ' num2str(dx)]});
    saveas(f, 'images/kineticDamping.png')
    
    semilogy(time(2:end), histPotEnergy);
    title({[distribution ' Potential Energy']; ['dt = ' num2str(DT) ', N = ' num2str(N, '%g') ', dx = ' num2str(dx)]});
    saveas(f, 'images/potentialDamping.png')
    
    semilogy(time(2:end), histTotEnergy);
    title({[distribution ' Total Energy']; ['dt = ' num2str(DT) ', N = ' num2str(N, '%g') ', dx = ' num2str(dx)]});
    saveas(f, 'images/totalDamping.png')
end

function phaseSpace(xp, vp, NG, Title)
    % Create a plot of the phaseSpace
    vp = [-1; vp; 1];
    xp = [min(xp); xp; max(xp)];
    dx = (max(xp) - min(xp))/NG;
    g1=floor(xp/dx-.5)+1;
    
    dv = (max(vp) - min(vp))/NG;
    g2 = floor(vp/dv-.5)+1;
    n = histcn([g1,g2], NG, NG);
    n1 = n';
    n1(size(n,1) + 1, size(n,2) + 1) = 0;
    
    % Wrap around
    n1(:, 1) = n1(:, 1) + n1(:, end);
    n1(:, end) = n1(:, 1);
    
    xb = linspace(min(xp),max(xp),size(n,1)+1);
    yb = linspace(min(vp),max(vp),size(n,1)+1);
    h = pcolor(xb,yb,n1);
    shading interp
    colormap(hot);
    xlabel('x','Fontsize',14);
    ylabel('v','Fontsize',14);
    title(Title,'Fontsize',12,'FontWeight','bold');
end

function [count edges mid loc] = histcn(X, varargin)
    % function [count edges mid loc] = histcn(X, edge1, edge2, ..., edgeN)
    %
    % Purpose: compute n-dimensional histogram
    %
    % INPUT
    %   - X: is (M x N) array, represents M data points in R^N
    %   - edgek: are the bin vectors on dimension k, k=1...N.
    %     If it is a scalar (Nk), the bins will be the linear subdivision of
    %     the data on the range [min(X(:,k)), max(X(:,k))] into Nk
    %     sub-intervals
    %     If it's empty, a default of 32 subdivions will be used
    %
    % OUTPUT
    %   - count: n-dimensional array count of X on the bins, i.e.,
    %         count(i1,i2,...,iN) = cardinal of X such that
    %                  edge1(i1) <= X(:,i1) < edge1(i1)+1 and
    %                       ...
    %                  edgeN(iN) <= X(:,iN) < edgeN(iN)+1
    %   - edges: (1 x N) cell, each provides the effective edges used in the
    %     respective dimension
    %   - mid: (1 x N) cell, provides the mid points of the cellpatch used in
    %     the respective dimension
    %   - loc: (M x N) array, index location of X in the bins. Points have out
    %     of range coordinates will have zero at the corresponding dimension.
    %
    % DATA ACCUMULATE SYNTAX:
    %   [ ... ] = histcn(..., 'AccumData', VAL);
    %   where VAL is M x 1 array. Each VAL(k) corresponds to position X(k,:)
    %   will be accumulated in the cell containing X. The accumulate result
    %   is returned in COUNT.
    %   NOTE: Calling without 'AccumData' is similar to having VAL = ones(M,1)
    %
    %   [ ... ] = histcn(..., 'AccumData', VAL, 'FUN', FUN);
    %     applies the function FUN to each subset of elements of VAL.  FUN is
    %     a function that accepts a column vector and returns
    %     a numeric, logical, or char scalar, or a scalar cell.  A has the same class
    %     as the values returned by FUN.  FUN is @SUM by default.  Specify FUN as []
    %     for the default behavior.
    %
    % Usage examples:
    %   M = 1e5;
    %   N = 3;
    %   X = randn(M,N);
    %   [N edges mid loc] = histcn(X);
    %   imagesc(mid{1:2},N(:,:,ceil(end/2)))
    %
    % % Compute the mean on rectangular patch from scattered data
    %   DataSize = 1e5;
    %   Lat = rand(1,DataSize)*180;
    %   Lon = rand(1,DataSize)*360;
    %   Data = randn(1,DataSize);
    %   lat_edge = 0:1:180;
    %   lon_edge = 0:1:360;
    %   meanData = histcn([Lat(:) Lon(:)], lat_edge, lon_edge, 'AccumData', Data, 'Fun', @mean);
    %
    % See also: HIST, ACCUMARRAY
    % 
    % Bruno Luong: <brunoluong@yahoo.com>
    % Last update: 25/August/2011

    if ndims(X)>2
        error('histcn: X requires to be an (M x N) array of M points in R^N');
    end
    DEFAULT_NBINS = 32;

    AccumData = [];
    Fun = {};

    % Looks for where optional parameters start
    % For now only 'AccumData' is valid
    split = find(cellfun('isclass', varargin, 'char'), 1, 'first');
    if ~isempty(split)
        for k = split:2:length(varargin)
            if strcmpi(varargin{k},'AccumData')
                AccumData = varargin{k+1}(:);
            elseif strcmpi(varargin{k},'Fun')
                Fun = varargin(k+1); % 1x1 cell
            end
        end
        varargin = varargin(1:split-1);
    end

    % Get the dimension
    nd = size(X,2);
    edges = varargin;
    if nd<length(edges)
        nd = length(edges); % wasting CPU time warranty
    else
        edges(end+1:nd) = {DEFAULT_NBINS};
    end

    % Allocation of array loc: index location of X in the bins
    loc = zeros(size(X));
    sz = zeros(1,nd);
    % Loop in the dimension
    for d=1:nd
        ed = edges{d};
        Xd = X(:,d);
        if isempty(ed)
            ed = DEFAULT_NBINS;
        end
        if isscalar(ed) % automatic linear subdivision
            ed = linspace(min(Xd),max(Xd),ed+1);
        end
        edges{d} = ed;
        % Call histc on this dimension
        [dummy loc(:,d)] = histc(Xd, ed, 1);
        % Use sz(d) = length(ed); to create consistent number of bins
        sz(d) = length(ed)-1;
    end % for-loop

    % Clean
    clear dummy

    % This is need for seldome points that hit the right border
    sz = max([sz; max(loc,[],1)]);

    % Compute the mid points
    mid = cellfun(@(e) 0.5*(e(1:end-1)+e(2:end)), edges, ...
                  'UniformOutput', false);

    % Count for points where all coordinates are falling in a corresponding
    % bins
    if nd==1
        sz = [sz 1]; % Matlab doesn't know what is one-dimensional array!
    end

    hasdata = all(loc>0, 2);
    if ~isempty(AccumData)
        count = accumarray(loc(hasdata,:), AccumData(hasdata), sz, Fun{:});
    else
        count = accumarray(loc(hasdata,:), 1, sz);
    end

    return

end % histcn

