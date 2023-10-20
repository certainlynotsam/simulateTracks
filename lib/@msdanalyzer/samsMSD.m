%After importing your tracks using samimport.m check the number of reported
%tracks 'Found X tracks in the file', this number you need to input into
%the first slot of the following command:

ma = msdanalyzer(1, 'nm', 'sec');
ma = ma.addAll(tracks); %turns tracks into a suitable format. 

%plot tracks like in 'samimport.m'
[hps, ha] = ma.plotTracks;
ma.labelPlotTracks

%Conduct MSD calculation
ma = ma.computeMSD;
figure
ma.plotMeanMSD(gca, true);
ma.fitMeanMSD