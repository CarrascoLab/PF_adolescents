function StartTiltStrengths=DPF_ExoAttnGetStartingContrast(Observer);
eval(sprintf('Folder=dir(''data/%s'');',Observer));
fileIwant=Folder(end).name;
eval(sprintf('load(''data/%s/%s'')',Observer,fileIwant));
StartTiltStrengths(1)=stimulus.stair{1}.strength(end);
StartTiltStrengths(2)=stimulus.stair{2}.strength(end);
end