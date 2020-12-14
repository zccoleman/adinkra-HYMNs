classdef adinkra
    %ADINKRA Creates an object representation of an adinkra
    
    properties
        ADJMatrix   %The adjacency matrix of the adinkra, with 
                    %different numbers indicating different colors and
                    %negatives indicating dashed edges
        nodes       %The total number of nodes
        open        %The number of open nodes
        closed      %The number of closed nodes
        colors      %The number of colors of the adinkra
        LMatrices   %A 3-dimensional array containing the L-matrices of 
                    %each color.  LMatrices(:,:,n) is the LMatrix of the 
                    %n'th color.
        hand        %0 means left-handed (i.e. even # of colors, Banchoff 
                    %matrix has values in the top LEFT), 1 means 
                    %right-handed (odd # of colors)
    end
    
    methods
        function obj = adinkra(matr)
            %ADINKRA Constructs an instance of an adinkra from an input 
            %adjacency matrix and assigns relevant properties
            obj.ADJMatrix=matr;
            obj.open=size(matr,1);
            obj.closed=size(matr,2);
            obj.nodes=obj.open+obj.closed;
            obj.colors=abs(max(max(matr)));
            LMatrices=zeros(obj.open,obj.closed,obj.colors);
            for i=1:obj.colors
                obj.LMatrices(:,:,i)=LColBlock(obj,i);
            end
            obj.hand=mod(obj.colors,2);
        end
        
        function lcol = LColBlock(a,color)
            %LCOLBLOCK returns the adjacency matrix of the color-induced
            %subgraph of the adinkra a.  This method is only called in the
            %constructor, then the matrices are saved in properties.
            h=a.open;
            w=a.closed;
            lcol=zeros(h,w);
            for i=1:h
                for j=1:w
                    if a.ADJMatrix(i,j)==color
                        lcol(i,j)=1;
                    elseif a.ADJMatrix(i,j)==color*-1
                        lcol(i,j)=-1;
                    end
                end
            end
        end
        
        
        function M = liftMatrix(a,m,w)
            %LIFTMATRIX creates a node-lifting operator with m as
            %the coefficient and a varaiable number of word parameters.
            %If multiple word parameters are specified, a node-lifting
            %operator is created for each one, and their product is
            %returned (order does not matter since all are diagonal).
            %The final operator has dimensions appropriate for the adinkra a.
            d=a.open;
            W=size(w,2);
            M = sym('noderaisingOp',[d,d]);
            M = M*0+eye(d);
            for i=1:W
                bin=w(i);
                for j=(d-1):(-1):0
                    if (bin>=2^j)
                        M(j+1,j+1)=M(j+1,j+1)*m;
                        bin=bin-2^j;
                    %else
                    %    M(i+1,i+1)=1;
                    end
                end
            end
        end
        
        
        function newL = LTilde(a,color,mb,mf,wb,wf)
            %LTILDE creates the tilded L-matrix as defined in eq.4.6 of the
            %adinkra a for the color, with m as the coefficient and w being
            %the word operator corresponding to the number of lifted nodes.
            L = a.LMatrices(:,:,color);
            %s=size(varargin,2);
            %wb=zeros(s/2);
            %wf=zeros(s/2);
            %for i=1:s
            %    if (mod(i,2)==1)
            %        wb(fix(i/s)+1)=varargin{i};
            %    else
            %        wf(i/2)=varargin{i};
            %    end
            %end
            newL = liftMatrix(a,mb,wb)*L*liftMatrix(a,1/mf,wf);
        end
        
        
        function newR = RTilde(a,color,mub,muf,wb,wf)
            %RTILDE creates a tilded R-matrix as defined in eq.4.6 for the
            %adinkra a, color, n as the coefficients, and w as the number
            %of raised nodes.  This matrix is equivalent to the transpose
            %of the tilded L-matrix with n = m^{-1}.
            R=a.LMatrices(:,:,color).';
            %s=size(varargin,2);
            %wb=zeros(s/2);
            %wf=zeros(s/2);

            newR=liftMatrix(a,muf,wf)*R*liftMatrix(a,1/mub,wb);

        end
        
        
        function cmatr = CBlock(a,color,mb,mub,mf,muf,wb,wf)
            %CBLOCK returns a "color matrix", which has zeros in the top
            %left and bottom right corners, and the L and R matrices in the
            %top right and bottom left corners respectively.  Input takes
            %the lifting operators Bm and Bn.  (Originally, input was just
            %the lifting coefficients m and n.  This change saves time by
            %not recalculating the operators with each call to CBlock.)
            d=a.nodes/2;
            cmatr = sym('cmatr',[d,d]);
            cmatr = cmatr*0;
            cmatr(1:d,(d+1):2*d)=LTilde(a,color,mb,mf,wb,wf);       %Bm*a.LMatrices(:,:,color);    
            cmatr((d+1):2*d,1:d)=RTilde(a,color,mub,muf,wb,wf);     %(Bn*a.LMatrices(:,:,color)).';
        end
        
        
        function bmatr = Banchoff(a,mb,mub,mf,muf,wb,wf)
            %BANCHOFF creates the Banchoff matrix for the adinkra using the
            %lifting coefficients m and n and the word parameter w.
            I=a.colors;
            bmatr=CBlock(a,I,mb,mub,mf,muf,wb,wf);
            for i=(I-1):-1:1
                bmatr=bmatr*CBlock(a,i,mb,mub,mf,muf,wb,wf);
            end
        end
        
        
        function eigs = HYMNs(a,varargin)
            %HYMNS makes a Banchoff matrix with the symbolic values m, n,
            %and r=m/n, then calculates the eigenvalues of the desired
            %parts of the matrix, depending on the type of adinkra.  Input
            %is the adinkra a followed by a list of word parameters.  There
            %must be an even number of word parameters.  They alternate
            %between bosonic and fermionic, starting from the first, which
            %should be bosonic.  Additional arguments of 0 do nothing (but
            %are sometimes necessary to ensure an even-sized list).
            mb=sym('mb');
            rhoB=sym('rhoB');
            mub=mb/rhoB;
            muf=sym('muf');
            rhoF=sym('rhoF');
            mf=muf/rhoF;
            s=size(varargin,2);
            wb=zeros(1,s/2);
            wf=zeros(1,s/2);
            for i=1:s
                if (mod(i,2)==1)
                    wb(fix(i/2)+1)=varargin{i};
                else
                    wf(i/2)=varargin{i};
                end
            end
            
            B=Banchoff(a,mb,mub,mf,muf,wb,wf);
            d=size(B,1)/2;
            if (a.hand==0) %want the eigenvalues of the top left and bottom right corners
                eigs=sym('eigs',[2*d,1]);
                eigs(1:d,1)=eig(B(1:d,1:d)); %Top left corner
                eigs(d+1:2*d,1)=eig(B(d+1:2*d,d+1:2*d)); %Bottom right corner
            else %want the eigenvalues of TR*BL and BL*TR
                eigs=sym('eigs',[d,2]);
                eigs(:,1)=eig(B(1:d,d+1:2*d)*B(d+1:2*d,1:d)); %TR*BL corners
                eigs(:,2)=eig(B(d+1:2*d,1:d)*B(1:d,d+1:2*d)); %BL*TR corners
                %eigs(a.open+1:a.nodes,1)=eig(BL(a,B));
            end
        end
        
        
        function eigs = DashlessHYMNs(a,wb,wb2,wf)
            %DASHLESSHYMNS creates an adinkra without dashings, then
            %calculates the HYMNs.
            b=adinkra(abs(a.ADJMatrix));
            eigs = HYMNs(b,wb,wb2,wf);
            %m=sym('m');
            %r=sym('r');
            %n=m/r;
            %B=Banchoff(b,m,n,w);
            %d=size(B,1)/2;
            %if (a.hand==0) %TL and BR
            %    eigs=sym('eigs',[2*d,1]);
            %    eigs(1:d,1)=eig(B(1:d,1:d));
            %    eigs(d+1:2*d,1)=eig(B(d+1:2*d,d+1:2*d));
            %else %TR and BL
            %    eigs=sym('eigs',[d,2]);
            %    eigs(:,1)=eig(B(1:d,d+1:2*d)*B(d+1:2*d,1:d)); %(1:a.open,1)
            %    eigs(:,2)=eig(B(d+1:2*d,1:d)*B(1:d,d+1:2*d));
                %eigs(a.open+1:a.nodes,1)=eig(BL(a,B));
            %end
        end
        
        
    end
end

