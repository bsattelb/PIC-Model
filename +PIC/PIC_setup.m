function PIC = PIC_setup(DT, NT, NG, N, distribution)
    % Returns a particle-in-cell evaluator that convertes a parameter array
    % to some quantity of interest.  In theh parameter array, the spatial
    % grid length, L, should be the first value.
    % DT is the temporal step
    % NT is the number of temporal steps
    % NG is the number of spatial gridpoints
    % N is the number of computational particles
    % distribution is the desired distribution
    %   This uses the files 'Initilization/initilization_<distribution>.m'
    %   and 'QOI_calc/calc_<distribution>.m' to initiliaze the particle
    %   positions and calculate the quantity of interest, respectively.
    %   Files for new distributions should be placed in those folders.  It
    %   is assumed that initilizations use the parameter vector and the 
    %   number of particles and the QOI calculation uses a vector of time 
    %   points and the 2nd mode of the spectral norm as parameters.
    
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
    
    % Create function handles for initilization and QOI calculation
    % Seems to need to be set up in this way to prevent reference issues -
    % simply using str2func and creating anonymous functions resulted in
    % issues
    function [xp, vp, rho_back, Q, QM] = initialization(params, N)
       curDir = pwd;
       cd([this_dir 'Initilization']);
       temp = ['initilization_' distribution];
       [xp, vp, rho_back, Q, QM] = feval(temp, params, N);
       cd(curDir);
    end
    
    function damprate = damp_calc(time, histEnergy)
        curDir = pwd;
        cd([this_dir 'QOI_calc']);
        temp = ['calc_' distribution];
        damprate = feval(temp, time, histEnergy);
        cd(curDir);
    end
    
    % Compute some unchanging values for use in the PIC computations
    p=1:N;
    p=[p p];
    un=ones(NG-1,1);
    Poisson=spdiags([un -2*un un],[-1 0 1],NG-1,NG-1);
    time = 0:DT:NT*DT;
    
    % Return a function handle to the PIC evaluator
    PIC = @fast_PIC;
    
    % Numerically calculate particle behavior
    function damprate = fast_PIC(params)
        L = params(1);
        dx = L/NG;
        [xp, vp, rho_back, Q, QM] = initialization(params, N);
        
        % history vectors and matrices
        histSpectrum = [];

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

           % take the Fourier transform of the electric field on the grid
           NFFT = 2^nextpow2(length(Eg)); % Next power of 2 from length of Eg
           Y = fft(Eg,NFFT)/length(Eg);
           histSpectrum = [histSpectrum 2*abs(Y(1:NFFT/2+1))];

        end

        % Find QOI using Spectral norm
        damprate = damp_calc(time, histSpectrum(2,:));
    end
end