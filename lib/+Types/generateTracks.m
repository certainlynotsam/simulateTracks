function tracks = generateTracks(motion_type, Diff, dT, N_particles, N_time_steps, N_dim, Size)

switch motion_type 

    case 'brownian'

    k = sqrt(2 * Diff * dT); % standard deviation from gaussian distribution
    
    tracks = cell(N_particles, 1); % store tracks in a cell array
    
    for i = 1 : N_particles
    
        % Time
        time = (0 : N_time_steps)' * dT;
    
        % Initial position
        X0 = [0 0 0]; %Size .* rand(1, N_dim);
    
        % Integrate uncorrelated displacement
        dX = k * randn(N_time_steps, N_dim);
        %dX(1, :) = X0;
        X = cumsum(dX, 1);
        X = cat(1,X0,X) + (rand(1, N_dim)*Size);
    
        % Store
        tracks{i} = [time X];
    
    end

    case 'confined' % random displacement + forced displacement

    kT = 4.2821e-21; % kBoltzman x T @ 37ºC

    k = sqrt(2 * Diff * dT);

    Ltrap = 0.100; % µm
    Ktrap = kT / Ltrap^2; % = thermal energy / trap size ^ 2

    tracks = cell(N_particles, 1);

    for i = 1 : N_particles % simulating over particles
    
        % Time
        time = (0 : N_time_steps-1)' * dT;
    
        % Initial position
        X0 = Size .* rand(1, N_dim);
    
        % Energy potential:
        V = @(x) 0.5 * Ktrap * sum (x .^ 2); % Unused, just to show
        Fx = @(x) - Ktrap * (x - X0); % Is a vector
    
        % Position
        X = zeros(N_time_steps, N_dim);
    
        % Init first step
        X(1, :) = X0;
    
        % Iterate
        for j = 2 : N_time_steps % simulating over points in a track
    
            dxtrap = Diff/kT * Fx(X(j-1,:)) * dT; % ad hoc confined disp
            dxbrownian = k * randn(1, N_dim); % ad hoc brownian disp
    
            X(j,:) = X(j-1,:) + dxtrap + dxbrownian; % add disp to prev track
    
        end
    
        % Store
        tracks{i} = [time X];
    
    end

    case 'mixed' % mixed brownian and confined motion

    fraction_brownian = 0.2; % proportion of tracks under brownian motion

    k = sqrt(2 * Diff * dT); % standard deviation from gaussian distribution

    kT = 4.2821e-21; % kBoltzman x T @ 37ºC

    Ltrap = 0.100; % µm

    Ktrap = kT / Ltrap^2; % = thermal energy / trap size ^ 2
    
    tracks = cell(N_particles, 1); % store tracks in a cell array
    
    for i = 1 : fraction_brownian*N_particles
    
        % Time
        time = (0 : N_time_steps)' * dT;
    
        % Initial position
        X0 = [0 0 0]; %Size .* rand(1, N_dim);
    
        % Integrate uncorrelated displacement
        dX = k * randn(N_time_steps, N_dim);
        %dX(1, :) = X0;
        X = cumsum(dX, 1);
        X = cat(1,X0,X) + (rand(1, N_dim)*Size);
    
        % Store
        tracks_brownian{i} = [time X];
    
    end

    for i = 1 : N_particles - numel(tracks_brownian) % simulating over particles
    
        % Time
        time = (0 : N_time_steps-1)' * dT;
    
        % Initial position
        X0 = Size .* rand(1, N_dim);
    
        % Energy potential:
        V = @(x) 0.5 * Ktrap * sum (x .^ 2); % Unused, just to show
        Fx = @(x) - Ktrap * (x - X0); % Is a vector
    
        % Position
        X = zeros(N_time_steps, N_dim);
    
        % Init first step
        X(1, :) = X0;
    
        % Iterate
        for j = 2 : N_time_steps % simulating over points in a track
    
            dxtrap = Diff/kT * Fx(X(j-1,:)) * dT; % ad hoc confined disp
            dxbrownian = k * randn(1, N_dim); % ad hoc brownian disp
    
            X(j,:) = X(j-1,:) + dxtrap + dxbrownian; % add disp to prev track
    
        end
    
        % Store
        tracks_confined{i} = [time X];
    end

    tracks = cat(2,tracks_brownian,tracks_confined);

end
