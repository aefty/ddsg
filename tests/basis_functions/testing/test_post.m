clc; close all;

a3 = sgParse;
a31 = sgParse;

interplateFile = '../interpolate.plt';
interplateFile1 = '../interpolate_post.plt';
originalFile ='../original.plt';
gridFile='../grid.plt';
surplusFile = '../surplus.plt';
basisFile ='../baisFunction.plt';

a3.loadData(...
    interplateFile,...
    originalFile,...
    gridFile,...
    surplusFile,...
    basisFile);


a31.loadData(...
    interplateFile1,...
    originalFile,...
    gridFile,...
    surplusFile,...
    basisFile);



disp('Total Error = Sum(Sum((Interplate_post - Interplate_pre).^2))'); 
sum(sum(a31.interpolant - a3.interpolant).^2)

disp('Average Error = mean(mean((Interplate_post - Interplate_pre).^2))');
mean(mean(a31.interpolant - a3.interpolant).^2)

disp('Max Error = max(max(((Interplate_post - Interplate_pre).^2))');
max(max(a31.interpolant - a3.interpolant).^2)

imagesc(a31.interpolant - a3.interpolant);                                        
colorbar


            
            




