% Block by Block Overall Accuracy comapred to Contrast
% First run DPF_ExoAttnAnalysis, then run this program


BlockAccuracyContrast = [];

totalBlocks = size(Blocks,1);
totalTrials = 48;

startTrial = 1;

for i = 1:totalBlocks    
    
    endTrial = 48*i;
    BlockAccuracyContrast(2,i) = (sum(temp3(startTrial:endTrial))/totalTrials);
    BlockAccuracyContrast(1,i) = Contrast(i);
    startTrial = startTrial + 48;
    
end

   
   