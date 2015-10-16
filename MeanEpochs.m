% Author: Rodolfo Abreu, ISR/IST, Universade de Lisboa, 2015

function mean_epochs = MeanEpochs(epocheddata)

N = size(epocheddata, 1);
L = size(epocheddata, 2);
mean_epochs = zeros(N, L);

for i = 1:N
    mean_epochs(i, :) = mean(squeeze(epocheddata(i, :, :)), 2);
end

return;
       
    
    
    

