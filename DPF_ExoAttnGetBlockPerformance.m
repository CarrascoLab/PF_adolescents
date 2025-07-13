function IndContrast=DPF_ExoAttnGetBlockPerformance(Observer,IndContrast);
eval(sprintf('Folder=dir(''data/%s'');',Observer));
fileIwant=Folder(end).name;
eval(sprintf('load(''data/%s/%s'')',Observer,fileIwant));
BlockData1=getTaskParameters(myscreen,task);
if BlockData1.nTrials==48
    if IndContrast==0
        IndContrast=stimulus.contrasts;
    else
        BlockData2=[BlockData1.randVars.targetOrientation' BlockData1.response'];
        BlockData3=BlockData2(:,1)==BlockData2(:,2);
        PropCorrect=mean(BlockData3);
        DiffCorrect=PropCorrect-.8;
        IndContrast=IndContrast-(IndContrast*DiffCorrect);
        if IndContrast>1
            IndContrast=1;
        end
    end
end