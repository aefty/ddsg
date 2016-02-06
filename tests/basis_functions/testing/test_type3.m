clc; close all;

a3 = sgParse;
    
interplateFile = 'type3/interpolate.plt';
originalFile ='type3/original.plt';
gridFile='type3/grid.plt';
surplusFile = 'type3/surplus.plt';
basisFile ='type3/baisFunction.plt';

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



            
            




