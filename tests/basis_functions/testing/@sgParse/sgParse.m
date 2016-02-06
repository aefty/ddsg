classdef sgParse < handle
    
    properties (Access = public)
        grid;
        interpolant;
        surplus;
        exact_solution;
        schematic;
        basis;
        
        dim;
        nno;
        df;
        level;
        gridType;

    end

    methods (Access = public)
        
        % Constructor
        function this = sgParse(this)
        end
          
        % Load Data
        function this = loadData(this,interplateFile,originalFile,gridFile,surplusFile,basisFile)

            this.interpolant = dlmread(interplateFile);
            this.exact_solution = dlmread(originalFile);
            
            temp = dlmread(gridFile,'\t',2);
            temp(:,3)=[];
            this.grid= temp;
            
            temp = [];
            
            temp = dlmread(surplusFile);
            this.dim = temp(1,1);
            this.nno = temp(1,2);
            this.df = temp(1,3);
            this.level = temp(1,4)+1;
            
           
            %Get gridType from surplus
            try
                this.gridType = temp(1,5);
            catch
                this.gridType = 1;
            end
            
            %Correct Type
            if(this.gridType==0)
                this.gridType=1;
            end
            
            temp(1,:) =[];            
            this.surplus =temp;
            
            % Load Bais Function
            data = dlmread(basisFile);
            this.basis = {};
            
            lvl = max(data(:,1));
            
            for i = 1:lvl
                this.basis{i} = [];
            end
            
            for row = data'
                this.basis{row(1)}(end+1,:) = [row(2),row(3)];
            end
            
           
        end

        % Visulaize Function & Error
        function this = visError(this)
            h = size(this.interpolant);

            [X,Y] = meshgrid(0:1/max(h):1-1/max(h));

            figure;
            subplot(3,2,1);
            a= flipud(this.interpolant);
            imagesc(a);
            title('Interploate');

            subplot(3,2,2);
            mesh(X,Y,this.interpolant);
            title('Interploate');


            subplot(3,2,3);
            a= flipud(this.exact_solution);
            imagesc(a);
            title('Original');

            subplot(3,2,4);
            mesh(X,Y,this.exact_solution);
            title('Original');
            
            
            delta = (this.interpolant-this.exact_solution);


            subplot(3,2,5);
            a= flipud(delta);
            imagesc(a);
            title('Error');

            subplot(3,2,6);
            mesh(X,Y,delta);
            title('Error');
            
            % Plot Error without Boundary

            figure
            delta(:,1)=[];
            delta(:,end)=[];
            delta(1,:)=[];
            delta(end,:)=[];
                        
            a= flipud(delta);
            imagesc(a);
            title('Error WO Boundary');  
            

            
                        
        end
        
        % Visulaize Schem
        function this = visSchem(this)
            this.getSchem;
            figure 
            for i = 1:this.level
                for j = 1:this.level-i+1
                    subplot(this.level,this.level,(i-1)*this.level+j) 
                    imagesc(this.schematic{i}{j}.grid);                    
                end 
            end            
        end

         % Visulaize Basis
        function this = visBasis(this)
            
            lvl =max(size(this.basis));
            
            figure 
            for i = 1:lvl
                subplot(lvl,1,i) 
                scatter(this.basis{i}(:,1),this.basis{i}(:,2));                    
            end
                
        end
        
        function this = visGrid(this)
            figure
            scatter(this.grid(:,1),this.grid(:,2));
            axis([0 1 0 1]) ;
        end
    end

    methods (Access = private)
        
        % Pars Scheme
        function this = getSchem(this)
            
            grp={};
           
            %1:Clenshaw-Curtis  or 3:Flipup
            if(this.gridType == 1 )
                
                s = 2^(this.level-1)+1;

                
                counter=zeros(this.level,this.level);

                for row = this.surplus'
                    grp{row(1)}{row(2)}.data = [];                    
                    grp{row(1)}{row(2)}.alpha = zeros(s,s);
                    grp{row(1)}{row(2)}.grid = zeros(s,s);
                    counter(row(1),row(2)) = counter(row(1),row(2))+1;
                end
                
                for row = this.surplus'
                    
                    coord=[];
                    
                    for d = 1:this.dim
                        level = row(d); 
                        index = row(d+this.dim);
                        
                        if (level ==1)
                            coord(d) = 0.5;
                        else
                            m = 2^(level-1)+1;
                            coord(d) = (index - 1) / (floor(m-1));
                        end
                    end
                    
                    a=min(floor(coord(1)*s)+1,s);
                    b=min(floor(coord(2)*s)+1,s);
                    
                    grp{row(1)}{row(2)}.data = [grp{row(1)}{row(2)}.data;coord,row(5)];
                    grp{row(1)}{row(2)}.alpha(a,b) = row(5);
                    grp{row(1)}{row(2)}.grid(a,b) = 1;
                    
                   
                end
                
                
            elseif(this.gridType == 3)
                
                error('Flipup  Not Done');
                
                          
            %2: Chebyshev-Gauss-Lobatto 
            elseif(this.gridType==2)
                
                
                error(' Chebyshev-Gauss-Lobatto  Not Done');
                
                
                                
            else
                error('Basis Function type Not detected');
            end
              
            this.schematic = grp;
       
        end
    
    end
end % class