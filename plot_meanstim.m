%%
figure
meanstim = mean(stim,2);
stemeanstim = abs(std(meanstim)/sqrt(size(meanstim,1)));
hist(meanstim)
biascheck = abs(mean(meanstim)) - abs(std(meanstim)/sqrt(size(meanstim,1)));
nsigma = biascheck/stemeanstim;
if biascheck<0
    title({['Mean of means of stimulus: ' num2str(mean(meanstim)) '+-' num2str(stemeanstim)], ['Bias: OK (biascheck=' num2str(biascheck) ')']})
else
    title({['Mean of means of stimulus: ' num2str(mean(meanstim)) '+-' num2str(stemeanstim)], ['Bias: NOT OK (biascheck=' num2str(biascheck) ';' num2str(nsigma) ' sigma)']})
end
xlabel('Mean of stimulus')
ylabel('Counts')
%%
figure
meannoisy = mean(stim_noisy,2);
stemeannoisy = abs(std(meannoisy)/sqrt(size(meannoisy,1)));
hist(meannoisy)
biascheck = abs(mean(meannoisy)) - stemeannoisy;
nsigma = biascheck/stemeanstim;
if biascheck<0
    title({['Mean of means of noisy stimulus: ' num2str(mean(meannoisy)) '+-' num2str(stemeannoisy)], ['Bias: OK (biascheck=' num2str(biascheck) ')']})
else
    title({['Mean of means of noisy stimulus: ' num2str(mean(meannoisy)) '+-' num2str(stemeannoisy)], ['Bias: NOT OK (biascheck=' num2str(biascheck) ';' num2str(nsigma) ' sigma)']})
end
xlabel('Mean of noisy stimulus')
ylabel('Counts')
%%
figure
meanstimuli = mean(stimuli,2);
stemeanstimuli = abs(std(meanstimuli)/sqrt(size(meanstimuli,1)));
hist(meanstimuli)
if mean(meanstimuli)>0.5
    biascheck = abs(mean(meanstimuli)) - stemeanstimuli;
    nsigma = (abs(mean(meanstimuli)) - 0.5)/stemeanstimuli;
else
    biascheck = abs(mean(meanstimuli)) + stemeanstimuli;
    nsigma = (0.5 - abs(mean(meanstimuli)))/stemeanstimuli;
end
if biascheck<0.5
    title({['Mean of means of belief: ' num2str(mean(meanstimuli)) '+-' num2str(stemeanstimuli)], ['Bias: OK (biascheck=' num2str(biascheck) ')']})
else
    title({['Mean of means of belief: ' num2str(mean(meanstimuli)) '+-' num2str(stemeanstimuli)], ['Bias: NOT OK (biascheck=' num2str(biascheck) ';' num2str(nsigma) ' sigma)']})
end
xlabel('Mean of belief')
ylabel('Counts')
%%
figure
meanstim2 = nanmean(stim2,2);
stemeanstim2 = abs(std(meanstim2)/sqrt(size(meanstim2,1)));
hist(meanstim2)
biascheck = abs(mean(meanstim2)) - abs(std(meanstim2)/sqrt(size(meanstim2,1)));
nsigma = biascheck/stemeanstim2;
if biascheck<0
    title({['Mean of means of experienced stimulus: ' num2str(mean(meanstim2)) '+-' num2str(stemeanstim2)], ['Bias: OK (biascheck=' num2str(biascheck) ')']})
else
    title({['Mean of means of experienced stimulus: ' num2str(mean(meanstim2)) '+-' num2str(stemeanstim2)], ['Bias: NOT OK (biascheck=' num2str(biascheck) ';' num2str(nsigma) ' sigma)']})
end
xlabel('Mean of experienced stimulus')
ylabel('Counts')