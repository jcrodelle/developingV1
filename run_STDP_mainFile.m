%% run_STDP_mainFile.m created by Jennifer Crodelle on April 22, 2020.
%% this code will run the c++ file title STDP_mainFile.cpp 
%% (which calls STDP_function.hpp) and save data files 
%% those data files will be read in and some weights plotted

clear all
N = 256;% number of cortical neurons -- a square number
Tfin = 100.0;  % in seconds
% gmax_cortical = 0.025;  %for rad connect
gmax_cortical = 0.003;  % for all-to-all
seed = 2000;
probGJs = 0.02;
cortLearnTime = 50;
flagForRadius = 0; %1 --> radius in cortex, 0 --> all-to-all

executable = 'STDP_mainFile';
command1 = ['rm ', executable];
system(command1)
command2 = 'g++ -lm -Wall -O2 STDP_mainFile.cpp -o STDP_mainFile';
system(command2);
path = '';

for j = 1:length(cortLearnTime)
    timeForSynapses = cortLearnTime(j);
    a = [N, Tfin*1000.0, cortLearnTime*1000.0, seed, gmax_cortical probGJs flagForRadius];
    B = ['./', num2str(executable),' ', num2str(a)]
    command3 = B;   % run code
    system(command3)
    
    endName = ['_gitCode'];
    name = [num2str(Tfin),'s',endName];
    movefile('LGNweights_final.csv',[path,'LGNweights_final_',name,'.csv']);
    movefile('weightMatrix_synapses.csv', [path,'weightMatrix_synapses_',name,'.csv']);
    movefile('weightMatrix_cortex.csv', [path,'weightMatrix_cortex_',name,'.csv']);
    movefile('parameters.txt', [path,'parameters_',name,'.txt']);
    movefile('neuronSpTimes.csv',[path,'neuronSpTimes_',name,'.csv']);
    movefile('electricConnections.csv',[path,'electricConnections_',name,'.csv']);
end

%% read in files and plot things

spikeTimes = dlmread([path,['neuronSpTimes_',name,'.csv']]);
neuronType = spikeTimes(:,1);
sisterID = spikeTimes(:,2);
numSpikes = spikeTimes(:,3);
spTimes = spikeTimes(:,4:end);
num_exc = length(find(neuronType==0));
I_inhib = find(neuronType==1);
I_exc = find(neuronType==0);

% load GJ connections:
eConn = dlmread([path,['electricConnections_',name,'.csv']]);
numberGJconn = eConn(:,1);
GJcolors = {(1/255)*[0 191 191],(1/255)*[120 120 120]};

% load LGN weights
WeightMatrix = dlmread([path,['weightMatrix_synapses_',name,'.csv']]);
[m,n] = size(WeightMatrix);
numMats = m/N;
M=sqrt(N);


% load cortical weights
WeightMatrix_cortex = dlmread([path,['weightMatrix_cortex_',name,'.csv']]); 

%% LGN weights
M1 = [];
figure
hold on
for i = 1:numMats-1
    Wmat = WeightMatrix(1+(i-1)*N:i*N,:);
    imagesc(Wmat'./0.02)
    title(['LGN weight: t = ', num2str(i*10), 'seconds']);
    colorbar
    caxis([0 1])
    axis([0.5 N+0.5 0.5 1000+0.5])
    M1 = [M1,getframe(gcf)];
    set(gca,'FontSize',15.0)
    xlabel('Receiving cell')
    ylabel('Sending cell')
end

%% cortical weights
exN = I_exc(40);
figure
hold on
aj = NaN(M,M);
M2 = [];
title('Cortical weights')
for j = 1:(numMats-1)
    Wmat2 = WeightMatrix_cortex(1+(j-1)*N:j*N,:)./gmax_cortical;
    newWmat = Wmat2(exN,:);
    newWmat(exN) = 2;
    for ii = 1:M
        aj(ii,:) = newWmat(1+(ii-1)*M:ii*M);
    end
    imagesc(aj)
    title(['Cortical weights: t = ', num2str(j*10), 'seconds']);
    colorbar
    caxis([0 2])
    M2 = [M2,getframe(gcf)];
    axis([0.5 M+0.5 0.5 M+0.5])
    set(gca,'FontSize',20.0)
    xlabel('Receiving cell')
    ylabel('Sending cell')
end

%% plot some spike times
figure
hold on
for jj = 1:N
    if neuronType(jj)==1
        scatter((1/1000).*spTimes(jj,1:numSpikes(jj)),jj*ones(numSpikes(jj),1),5,'b')
        hold on
    else
        scatter((1/1000).*spTimes(jj,1:numSpikes(jj)),jj*ones(numSpikes(jj),1),5,'r')
        hold on
    end
end
axis([580 581 0 N])
set(gca,'FontSize',20.0)
title('Raster plot')
xlabel('Time (ms)')
ylabel('Neuron index')

