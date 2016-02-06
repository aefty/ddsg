clc; close all;

a3 = sgParse;

interplateFile = '../interpolate.data';
originalFile ='../exact.data';
gridFile='../grid.data';
surplusFile = '../surplus.data';
basisFile ='../baisFunction.data';

a3.loadData(...
    interplateFile,...
    originalFile,...
    gridFile,...
    surplusFile,...
    basisFile);


a3.visError();
%a3.visSchem();
a3.visBasis();
a3.visGrid();



            
            




