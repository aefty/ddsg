clc; close all;

a1 = sgParse;
    
interplateFile = 'type1/interpolate.plt';
originalFile ='type1/original.plt';
gridFile='type1/grid.plt';
surplusFile = 'type1/surplus.plt';
basisFile ='type1/baisFunction.plt';

a1.loadData(...
    interplateFile,...
    originalFile,...
    gridFile,...
    surplusFile,...
    basisFile);

a1.visError();
a1.visSchem();
a1.visBasis();
a1.visGrid();


            
            




