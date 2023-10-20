function plotTracks(tracks, N_particles,space_units,n)

colors = turbo(N_particles);

figure
hold('on');
hps = NaN(N_particles, 1);

for i = 1 : N_particles
        
    track = tracks{i};
    trackName = sprintf('Track %d', i );
    
    x = track(:,2);
    y = track(:,3);
    z = track(:,4);
    
    hps(i) =  plot3(x, y, z, 'Color', colors(i,:), 'DisplayName', trackName ,'linewidth',1.5);
        
end

xlabel(['x (' space_units ')']); 
ylabel(['y (' space_units ')']); 
zlabel(['z (' space_units ')']);
view([45 45 45])