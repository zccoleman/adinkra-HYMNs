classdef adinkra
    %ADINKRA Creates an object representation of an adinkra and performs
    %necessary calculations on the adinkra.
    %%%%%%%%%%%%   Detailed explanation goes here
    
    properties
        ADJMatrix    %The Vertex-edge incidence matrix of the adinkra, with 
                    %different numbers indicating different colors
        nodes      %The total number of nodes
        open       %The number of open nodes
        closed     %The number of closed nodes
        colors     %The number of colors of the adinkra
        LMatrices  %A 3-dimensional array containing the L-matrices of each color.  LMatrices(:,:,n) is the LMatrix of the n'th color.
        hand %0 means left-handed (i.e. even # of colors, Banchoff matrix has values in the top LEFT), 1 means right-handed (odd # of colors)
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
            %the coefficient and w as the word parameter describing the number of lifted nodes.
            %The operator has size appropriate for the adinkra a.
            d=a.open;
            M = sym('noderaisingOp',[d,d]);
            M = M*0;
            bin=w;
            for i=d-1:-1:0
                if (bin>=2^i)
                    M(i+1,i+1)=m;
                    bin=bin-2^i;
                else
                    M(i+1,i+1)=1;
                end
            end
        end
        
        function newL = LTilde(a,color,M) %Defunct
            %LTILDE creates the tilded L-matrix as defined in eq.4.6 of the
            %adinkra a for the color, with m as the coefficient and w being
            %the word operator corresponding to the number of lifted nodes.
            L = a.LMatrices(:,:,color);
            newL = M*L;
        end
        
        function newR = RTilde(a,color,M) %Defunct
            %RTILDE creates a tilded R-matrix as defined in eq.4.6 for the
            %adinkra a, color, n as the coefficients, and w as the number
            %of raised nodes.  This matrix is equivalent to the transpose
            %of the tilded L-matrix with n = m^{-1}.
            newR=LTilde(a,color,M).';

        end
        
        function cmatr = CBlock(a,color,Bm,Bn)
            %CBLOCK returns a "color matrix", which has zeros in the top
            %left and bottom right corners, and the L and R matrices in the
            %top right and bottom left corners respectively.  Input takes
            %the lifting operators Bm and Bn.  (Originally, input was just
            %the lifting coefficients m and n.  This change saves time by
            %not recalculating the operators with each call to CBlock.)
            d=a.nodes/2;
            cmatr = sym('cmatr',[d,d]);
            cmatr = cmatr*0;
            cmatr(1:d,(d+1):2*d)=Bm*a.LMatrices(:,:,color);    %LTilde(a,color,Bm);
            cmatr((d+1):2*d,1:d)=(Bn*a.LMatrices(:,:,color)).';  %RTilde(a,color,Bn);
        end
        
        function bmatr = Banchoff(a,m,n,w)
            %BANCHOFF creates the Banchoff matrix for the adinkra using the
            %lifting coefficients m and n and the word parameter w.
            I=a.colors;
            Bm=liftMatrix(a,m,w);
            Bn=liftMatrix(a,1/n,w);
            bmatr=CBlock(a,I,Bm,Bn);
            for i=(I-1):-1:1
                bmatr=bmatr*CBlock(a,i,Bm,Bn);
            end
        end
        
        
        function eigs = HYMNs(a,w)
            %HYMNS makes a Banchoff matrix with the symbolic values m, n,
            %and r=m/n, then calculates the eigenvalues of the desired part
            %of the matrix, according to the type of adinkra.
            m=sym('m');
            r=sym('r');
            n=m/r;
            B=Banchoff(a,m,n,w);
            d=size(B,1)/2;
            if (a.hand==0) %Then we want the eigenvalues of the top left and bottom right corners
                eigs=sym('eigs',[2*d,1]);
                eigs(1:d,1)=eig(B(1:d,1:d)); %Top left corner
                eigs(d+1:2*d,1)=eig(B(d+1:2*d,d+1:2*d)); %Bottom right corner
            else %Then we want the eigenvalues of TR*BL and BL*TR
                eigs=sym('eigs',[d,2]);
                eigs(:,1)=eig(B(1:d,d+1:2*d)*B(d+1:2*d,1:d)); %TR*BL corners
                eigs(:,2)=eig(B(d+1:2*d,1:d)*B(1:d,d+1:2*d)); %BL*TR corners
                %eigs(a.open+1:a.nodes,1)=eig(BL(a,B));
            end
        end
        
        function eigs = DashlessHYMNs(a,w)
            %DASHLESSHYMNS creates an adinkra without dashings, then
            %calculates the HYMNs.
            m=sym('m');
            r=sym('r');
            n=m/r;
            b=adinkra(abs(a.ADJMatrix));
            B=Banchoff(b,m,n,w);
            d=size(B,1)/2;
            if (a.hand==0) %TL and BR
                eigs=sym('eigs',[2*d,1]);
                eigs(1:d,1)=eig(B(1:d,1:d));
                eigs(d+1:2*d,1)=eig(B(d+1:2*d,d+1:2*d));
            else %TR and BL
                eigs=sym('eigs',[d,2]);
                eigs(:,1)=eig(B(1:d,d+1:2*d)*B(d+1:2*d,1:d)); %(1:a.open,1)
                eigs(:,2)=eig(B(d+1:2*d,1:d)*B(1:d,d+1:2*d));
                %eigs(a.open+1:a.nodes,1)=eig(BL(a,B));
            end
        end
    end
end

