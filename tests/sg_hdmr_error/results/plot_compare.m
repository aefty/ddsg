close all; clear;
dirFolder='';
tit = {'F1-OSCILLATORY','F2-PRODUCT PEAK','F3-CORNER PEAK','F4-LOG GAUSSIAN','F5-CONTINUOUS','F6-MAX'};

folderStruc=dir('test*');
folderName={};

FN = max(size(folderStruc));

for i = 1:FN
    folderName{i} = folderStruc(i).name;
end

value.exact=[];
value.sg=[];
value.hdmr=[];

points.sg=[];
points.hdmr=[];

err.sg=[];
err.hdmr=[];

%Funciton values -F
value.exact = csvread(strcat(dirFolder,folderName{i},'/exact.csv'));
value.sg    = csvread(strcat(dirFolder,folderName{i},'/sgInterpolate.csv'));
value.hdmr  = csvread(strcat(dirFolder,folderName{i},'/hdmrInterpolate.csv'));

%Input Value   X
X = csvread(strcat(dirFolder,folderName{i},'/x.csv'));

%Number of points in SG & HDMR
points.sg   = csvread(strcat(dirFolder,folderName{i},'/gridpoints_sg.csv'));
points.hdmr = csvread(strcat(dirFolder,folderName{i},'/gridpoints_hdmr.csv'));

%Error
err.sg   = value.sg - value.exact;
err.hdmr = value.hdmr - value.exact;

%Number of functions
N=size(value.exact,2);

%Dimention
D=size(X,2);

%Number of data points
ND=max(size(X));



for i=1:N
 RowName{i}=tit{i};%sprintf('F(%d)',i);
end

format short e
L1sg = sum(abs(err.sg),1)/ND;
L2sg = sqrt(sum(err.sg.^2,1))/ND;
LIsg = max(abs(err.sg));

L1hdmr = sum(abs(err.hdmr),1)/ND;
L2hdmr = sqrt(sum(err.hdmr.^2,1))/ND;
LIhdmr = max(abs(err.hdmr));


L1_error=[L1sg',L1hdmr',L1hdmr'-L1sg'];
L2_error=[L2sg',L2hdmr',L2hdmr'-L2sg'];
LI_error=[LIsg',LIhdmr',LIhdmr'-LIsg'];

disp('Error Tables [SG,HDMR,HDMR/SG-1]');


T = table(L1_error,L2_error,LI_error,...
    'RowNames',RowName)


%Genearte combination of dimentions in 2d
dimCuts = nchoosek((1:D),2);
Dset = size(dimCuts,1);

h=100;
[xq,yq] = meshgrid(0:1/h:1, 0:1/h:1);

disp('Generating Plots Interpolation ...');
for i=1:N
    figure
    plotCell=0;
    
    for j=1:Dset
        
        D1=dimCuts(j,1);
        D2=dimCuts(j,2);
        
        x=X(:,D1);
        y=X(:,D2);

        v0= value.exact(:,i);
        v1= value.sg(:,i);
        v2= value.hdmr(:,i);

        % Exact Solution           
        plotCell=plotCell+1;        
        subplot(Dset,3, plotCell);
        vq0 = griddata(x,y,v0,xq,yq);
        mesh(xq,yq,vq0);
        title(sprintf('Exact Solution-(D%d,D%d)',D1,D2));
        xlabel(sprintf('D%d',D1));
        ylabel(sprintf('D%d',D2));
        zlabel('Value');
        
        % SG Solution
        plotCell=plotCell+1;        
        subplot(Dset,3, plotCell);
        vq1 = griddata(x,y,v1,xq,yq);
        mesh(xq,yq,vq1);
        title(sprintf('SG Solution-(D%d,D%d)',D1,D2));
        xlabel(sprintf('D%d',D1));
        ylabel(sprintf('D%d',D2));
        zlabel('Value');
        
        % HDMR Solution                
        plotCell=plotCell+1;        
        subplot(Dset,3, plotCell);
        vq2 = griddata(x,y,v2,xq,yq);
        mesh(xq,yq,vq2);
        title(sprintf('HDMR Solution-(D%d,D%d)',D1,D2));
        xlabel(sprintf('D%d',D1));
        ylabel(sprintf('D%d',D2));
        zlabel('Value');
        
    end
     suplabel(sprintf('Solution %s',tit{i}),'t');
end

pause
close all


disp('Generating Plots of ERRORS ...');
for i=1:N
    figure
    plotCell=0;
    
    for j=1:Dset
        
        D1=dimCuts(j,1);
        D2=dimCuts(j,2);
        
        x=X(:,D1);
        y=X(:,D2);

        v1= abs(err.sg(:,i));
        v2= abs(err.hdmr(:,i));
        v3= abs(v2-v1);
        
        % SG ABS Error
        plotCell=plotCell+1;        
        subplot(Dset,3, plotCell);
        vq1 = griddata(x,y,v1,xq,yq);
        contourf(xq,yq,vq1);
        title(sprintf('SG ABS Error-(D%d,D%d)',D1,D2));
        xlabel(sprintf('D%d',D1));
        ylabel(sprintf('D%d',D2));
        zlabel('Value');
        colorbar
        
        % HDMR ABS Error                
        plotCell=plotCell+1;        
        subplot(Dset,3, plotCell);
        vq2 = griddata(x,y,v2,xq,yq);
        contourf(xq,yq,vq2); 
        title(sprintf('HDMR ABS Error-(D%d,D%d)',D1,D2));
        xlabel(sprintf('D%d',D1));
        ylabel(sprintf('D%d',D2));
        zlabel('Value');
        colorbar
        
        % HDMR error vs SG Error          
        plotCell=plotCell+1;        
        subplot(Dset,3, plotCell);
        vq3 = griddata(x,y,v3,xq,yq);
        contourf(xq,yq,vq3);    
        title(sprintf('HDMR error vs SG Error-(D%d,D%d)',D1,D2));
        xlabel(sprintf('D%d',D1));
        ylabel(sprintf('D%d',D2));
        zlabel('Value');
        colorbar
        
    end
     suplabel(sprintf('Error Function %s',tit{i}),'t');
end
