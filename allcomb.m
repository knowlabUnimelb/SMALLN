function A = allcomb(varargin) ;
% ALLCOMB - All combinations
%    B = ALLCOMB(A1,A2,A3,...,AN) returns all combinations of the elements
%    in A1, A2, ..., and AN. B is P-by-N matrix is which P is
%    prod(numel(A1),..,numel(AN)). Empty inputs are ignored.
%
%    Example:
%    allcomb([1 3 5],[-3 8],[],[0 1]) ;
%       1  -3   0
%       1  -3   1
%       1   8   0
%       ...
%       5  -3   1
%       5   8   0
%       5   8   1
%    
%    ALLCOMB(A1,..AN,'matlab') causes the first column to change fastest.
%    This is more consistent with matlab indexing. Example:
%    allcomb(1:2,3:4,5:6,'matlab') %-> 
%      1   3   5
%      2   3   5
%      1   4   5
%      ...
%      2   4   6
%
%    This functionality is also known as the cartesion product.
%
%    See also PERMS, 
%    and COMBN (Matlab Central FEX)

% for Matlab R13
% version 1.0 (feb 2006)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History
% 1.1 (feb 2006), removed minor bug in empty cell arrays; 
%     added option to let the first input run fastest (suggestion by JD)


error(nargchk(1,Inf,nargin)) ;

% check for empty inputs
q = ~cellfun('isempty',varargin) ;
if any(~q), warning('Empty inputs are ignored') ; end

ni = sum(q) ;

argn = varargin{end} ;
if strcmpi(argn,'matlab') || strcmpi(argn,'john'),
    % based on a suggestion by JD on the FEX
    ni = ni-1 ;
    ii = 1:ni ;
    q(end) = 0 ;
else
    % enter arguments backwards, so last one (AN) is changing fastest
    ii = ni:-1:1 ;
end

if ni==0,
    A = [] ;
else
    args = varargin(q) ;
    if ~all(cellfun('isclass',args,'double')),
        error('All arguments should be arrays of doubles') ;
    end
    if ni==1,
        A = args{1}(:) ;
    else
        
        [A{1:ni}] = ndgrid(args{ii}) ;
        % concatenate
        A = reshape(cat(ni+1,A{:}),[],ni) ;
        % flip if AN is changing fastest
        A = A(:,ii) ;
    end
end


