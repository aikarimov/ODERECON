function sigma = deglexord(varargin)
%DEGLEXORD returns degree-lexicographic order for monomials in a form of
%matrix of powers by each monomial term
%   Syntax:
%
%   sigma = DEGLEXORD(dmax,M)
%
%   sigma = DEGLEXORD(dmin, dmax,M)
%
%DEGLEXORD returns degree-lexicographic ordering sigma of order d and
%dimension M
%   Syntax:
%
%   dmin - min order <= 0; by default, dmin = 0
%   dmax - max order
%   M - dimension
%   sigma is [D x M] matrix
%   D - number of unique monomials, depends on the min and max order
%   M - dimension

%   example:
%
%   sigma = deglexord(3,3) returns
%   sigma = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 2 0 0; 1 1 0; 1 0 1; 0 2 0; 0 1 1; 0 0 2; 3 0 0; 2 1 0;  1 2 0;  2 0 1; 1 0 2; 0 3 0; 0 2 1; 0 1 2; 0 0 3];

if nargin == 2
    dmin = 0;
    dmax = varargin{1,1};
    M =    varargin{1,2};
end

if nargin == 3
    dmin = varargin{1,1};
    dmax = varargin{1,2};
    M =    varargin{1,3};
end
%S = zeros(1,M); % order 0
%for d = 1:dmax  % ordering from 1 to d

S = zeros(0,M); %empty 0 x M

Dneg = max([0,-dmin]);
Dpos = max([0,dmin]);

while (Dpos <= dmax) && (Dneg >= 0)% ordering from 1 to dmax - dmin
    dcur = Dneg + Dpos; %shifting param
    Ind = ones(dcur,1);
    
    %create matrix V, height equals to dcur, width is dimension M
    V = zeros(dcur,M);
    %make a column of -1 and 1 with respect to degree
    colsign = ones(dcur,1);
    colsign(1:Dneg) = -1;
    %add a columt as a first column of matrix V
    V(:,1) = colsign;
    
    
    while( abs(sum(abs(V(:,M)))) ~= dcur ) %while the last column is not full
        S = adds(V,S);
        V = zeros(dcur,M); %zeroing V once again    
        %shift indices
        Ind = shiftworward(Ind, dcur, M);
        for j = 1:dcur
            V(j,Ind(j)) = colsign(j);
        end
    end
    S = adds(V,S);
    %reduce degree of neg col
    if Dpos == dmax %reduce Dneg if Dpos is maimal, until the loop stops
        Dneg = Dneg - 1;
    end
    %increase degree of positive col
    if Dpos < dmax
        Dpos = Dpos + 1;
    end   
end
sigma = S;
end

function S = adds(V,S)
Vnew = sum(V,1);
if ~israwcontained(Vnew, S)
   S = [S; Vnew];
end
end

function ind = shiftworward(ind, dit, M)
if ind(dit) < M
    ind(dit) = ind(dit) + 1;
else
    if dit > 1
        ind = shiftworward(ind, dit - 1, M);
        ind(dit) = 1;
    end
end
end