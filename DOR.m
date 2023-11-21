function A_new = DOR(D, method)
%  ISPP_DOR, Decending order recovery
%  Input D: demand matrix, which is also a distance matrix here
%        Method: The method to build the graph
%                'classic' build a spanning tree first and then 
%                          remove redundant link directly from the demand
%                          complete graph
%                'advanced' remove redundant link directly from the demand
%                complete graph
%  Output: A_new : A weigthed adjacency matrix
    
if nargin < 1
    error('Not enough input arguments.');
elseif nargin > 2
    error('Too many input arguments.');
end
if nargin == 1
    %-----------------------
    method = 'advanced';
end

switch method
    case 'classic'
        G_D = graph(D,'lower');

        T = minspantree(G_D);
        A_T = full(adjacency(T,"weighted"));
        
        max_dif_distance_value =  max(max(tril(distances(T))-tril(D)));
        [maxrow,maxcol] = find(tril(distances(T)-D) == max_dif_distance_value,1);
        A_new = A_T;
        while max_dif_distance_value>0.00000001
       
            A_new(maxrow,maxcol) = D(maxrow,maxcol);
            A_new(maxcol,maxrow) = D(maxrow,maxcol);
            G_new = graph(A_new,'lower');
            A_new = full(adjacency(G_new,"weighted"));
            max_dif_distance_value =  max(max(distances(G_new)-D));
            [maxrow,maxcol] = find((distances(G_new)-D) == max_dif_distance_value,1);
        end    
        G1 = graph(A_new);
        S1 = distances(G1);
        A_new = remove_triangle(A_new,S1);
        
    case 'advanced'
        A_new = D;
        N = size(A_new,1);
        S_med = D;
        for i=1:N
            for j =i+1:N
                if A_new(i,j)~=0
                    for k = 1:N
                        if  k~=i && k~=j && A_new(i,j)>= S_med(i,k)+S_med(k,j)
                            A_new(i,j) =0;
                            A_new(j,i)=0;
                            break
                        end
                    end
                end
            end
        end
end
end


function A_OUT = remove_triangle(A_new,S1)
    
    N1 = size(A_new,1);
    % 
    S_med1 = S1;
    for i1=1:N1
        for j1 =i1+1:N1
            if A_new(i1,j1)~=0
                for k1 = 1:N1
                    if  k1~=i1 && k1~=j1 && A_new(i1,j1)>= S_med1(i1,k1)+S_med1(k1,j1)
                        A_new(i1,j1) =0;
                        A_new(j1,i1)=0;
                        break
                    end
                end
            end
        end
    end
    A_OUT =A_new; 
end








