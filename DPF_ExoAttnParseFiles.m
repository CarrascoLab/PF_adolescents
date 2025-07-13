function [TestBlocks]=DPF_ExoAttnParseFiles(Observer);
eval(sprintf('Folder=dir(''data/%s'');',Observer));
if size(Folder,1)
FileNum=length(Folder);
for File=3:FileNum
    DataIndex(File)=strcmp(Folder(File).name(end-2:end),'mat');
end
Counter=0;
for File=3:FileNum
    if DataIndex(File) && ~strcmp(Folder(File).name(1),'.')
        Counter=Counter+1;
        eval(sprintf('load(''data/%s/%s'')',Observer,Folder(File).name));
        if stimulus.trialend==48
            FullBlockIndex(Counter)=1;
        else
            FullBlockIndex(Counter)=0;
        end
        FileType{Counter,1}=task{1}.taskFilename;
        FileType{Counter,2}=Folder(File).name(1:6);
        FileType{Counter,3}=Folder(File).name(12:13);
        BlockPrep(Counter,1)=str2num(Folder(File).name(1:6));
        BlockPrep(Counter,2)=str2num(Folder(File).name(12:13));
    end
end
for i=1:length(FileType)
    if FullBlockIndex(i)
        checkFunction=strcmp(FileType{i,1},'DPF_ExoAttnExp.m')||strcmp(FileType{i,1},'DPF_ExoAttnExp_ExtendedTiming.m');
    TestIndex(i)=checkFunction;
    end
end
TestBlocks=BlockPrep(TestIndex,:);
else
    fprintf('Participant not found. Check Observer initials.\r')
TestBlocks=[];
end

