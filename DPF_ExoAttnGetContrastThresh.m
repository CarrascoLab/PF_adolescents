function [ReversalThresh]=DPF_ExoAttnGetContrastThresh(Observer);
eval(sprintf('Folder=dir(''data/%s'');',Observer));
TiltMat=1;
FileNum=length(Folder);
for File=3:FileNum
    DataIndex(File)=strcmp(Folder(File).name(end-2:end),'mat');
end
Counter=0;
for File=3:FileNum
    if DataIndex(File) && ~strcmp(Folder(File).name(1),'.')
        Counter=Counter+1;
        eval(sprintf('load(''data/%s/%s'')',Observer,Folder(File).name));
        FileType{Counter,1}=task{1}.taskFilename;
        FileType{Counter,2}=Folder(File).name(1:6);
        FileType{Counter,3}=Folder(File).name(12:13);
        BlockPrep(Counter,1)=str2num(Folder(File).name(1:6));
        BlockPrep(Counter,2)=str2num(Folder(File).name(12:13));
        if stimulus.trialend==48 
            if strcmp(FileType{Counter,1},'DPF_ExoAttnContrastThresh.m')
            fileIwant=Folder(File).name;
            end
        end
    end
end

eval(sprintf('load(''data/%s/%s'')',Observer,fileIwant));
Data{TiltMat}=getTaskParameters(myscreen,task);
StairCaseIndex{TiltMat}=Data{TiltMat}.randVars.staircase;
SC1counter=1;
SC2counter=1;
for i=1:length(StairCaseIndex{TiltMat})
    if StairCaseIndex{TiltMat}(i)==1
        StairCaseStength{TiltMat}(i)=stimulus.stair{1}.strength(SC1counter);
        SC1counter=SC1counter+1;
    else
        StairCaseStength{TiltMat}(i)=stimulus.stair{2}.strength(SC2counter);
        SC2counter=SC2counter+1;
    end
end
    ReversalThresh=(stimulus.stair{1}.threshold+stimulus.stair{2}.threshold)/2;
