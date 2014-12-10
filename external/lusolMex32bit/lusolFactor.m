function [L,U,p,q,options,stats] = lusolFactor(A,options)
%
%        options   = lusolSet;
%        [L,U,p,q,options] = lusolFactor(A,options);
%
% lusolFactor computes the sparse factorization A = L*U
% for a square or rectangular matrix A.  The vectors p, q
% are row and column permutations giving the pivot order.
%    L(p,p) is unit lower triangular with subdiagonals bounded
%           by options.FactorTol,
%    U(p,q) is upper triangular or upper trapezoidal,
%           depending on the shape and rank of A.
%    A(p,q) = L(p,p)*U(p,q) would be a truly triangular
%             (or trapezoidal) factorization.
%    The rank of A tends to be reflected in the rank of U,
%    especially if options.Pivoting = 'TRP' or 'TCP'.
%
% Example use:
%    options = lusolSet(options,'pivot','TRP','factor',2.0);
% or equivalently
%    options = lusolSet;
%    options.Pivoting   = 'TRP';
%    options.FactorTol  = 2.0;
%
%   [L,U,p,q,options] = lusolFactor(A,options);
%
%   inform  = options.Inform;
%   rank    = options.Rank;

% The F-Mex Matlab interface to LUSOL is maintained by
% Michael O'Sullivan and Michael Saunders, SOL, Stanford University.
%
% Known Bugs:
%
% options.KeepLU = 'No' does not work.
% The Mex file still tries to create L and U.
%
% 02 Feb 1999: MJO: Developed F-MEX interface to LUSOL's lu1fac.
% 18 Oct 2000: MAS: Added options structure.
% 15 Apr 2001: MJO: options is now an output parameter.
% 14 Aug 2002: MAS: Added TRP option.
% 29 Apr 2004: MAS: L and U are returned as given by lu1fac.
%              There is no need to permute them to triangular form
%              because Matlab handles permuted triangles correctly
%              in statements like x = U\(L\b);
% 07 Jul 2005: MJO: F-Mex lu1fac now works correctly on rectangular A.
% 18 Jan 2010: RMF: Interfaced to Nick Hendersons object oriented
%                   matlab interface to 64 bit lusol.

[m,n] = size(A);
if ~issparse(A), A = sparse(A); end     % Make sure A is sparse.
stats=[];
archstr = computer('arch');
archstr = lower(archstr);
switch archstr
    case 'glnx86'
        luparm    = zeros(30,1);                % Store options in LUSOL arrays.
        parmlu    = zeros(30,1);
        
        luparm(1) = options.PrintFile;
        luparm(2) = options.PrintLevel;
        luparm(3) = options.MaxCol;
        if strcmpi(options.Pivoting,'TPP'), luparm(6) = 0; end
        if strcmpi(options.Pivoting,'TRP'), luparm(6) = 1; end
        if strcmpi(options.Pivoting,'TCP'), luparm(6) = 2; end
        if strcmpi(options.KeepLU  ,'No' ), luparm(8) = 0; end
        if strcmpi(options.KeepLU  ,'Yes'), luparm(8) = 1; end
        
        parmlu(1) = options.FactorTol;
        parmlu(2) = options.UpdateTol;
        parmlu(3) = options.DropTol;
        parmlu(4) = options.Utol1;
        parmlu(5) = options.Utol2;
        parmlu(6) = options.Uspace;
        parmlu(7) = options.Dense1;
        parmlu(8) = options.Dense2;
        [L,U,p,q,luparm,parmlu] = lu1fac(A,luparm,parmlu);   % Factorize A = L*U
        
        % 29 Apr 2004: No need to do this:
        % if options.KeepLU
        %        L = L(p,p);                  % Make L and U strictly triangular
        %        U = U(p,q);
        % end
        
        inform         = luparm(10);
        options.Inform = inform;
        options.Rank   = luparm(16);
        options.Nsing  = luparm(11);
        options.Growth = parmlu(16);
        
        if inform > 0
            disp(['Note: lu1fac returned Inform = ' num2str(inform)])
        end
    case 'glnxa64'
        if exist('options','var')
                       
            % generate default options from: options = lusol.luset();
            %initialise
            mylu = lusol(1);
            options64 = mylu.luset();
            
%             The parameter map between lu1fac and lusol_mex is:
%             | LUSOL     | lu1fac    | lusol_mex |
%             | parmlu(1) | FactorTol | Ltol1     |
%             | parmlu(2) | UpdateTol | Ltol2     |
%             | parmlu(3) | DropTol   | small     |
%             | parmlu(4) | Utol1     | Utol1     |
%             | parmlu(5) | Utol2     | Utol2     |
%             | parmlu(6) | Uspace    | Uspace    |
%             | parmlu(7) | Dense1    | dens1     |
%             | parmlu(8) | Dense2    | dens2     |
            if isfield(options,'PrintFile')
                % options64.       = options.PrintFile;%?
            end
            if isfield(options,'PrintLevel')
                PrintLevel=options.PrintLevel;
                % options64.        = options.PrintLevel;%?
            else
                PrintLevel=0;
            end
            if isfield(options,'MaxCol')
                options64.maxcol = options.MaxCol;
            end
            if isfield(options,'Pivoting')
                options64.pivot = options.Pivoting;
            end
            if isfield(options,'KeepLU')
                if strcmp(options.KeepLU,'Yes')
                    options64.keepLU = 1;
                else
                    options64.keepLU = 0;
                end
            end
            if isfield(options,'FactorTol')
                options64.Ltol1 = options.FactorTol;
            end
            if isfield(options,'UpdateTol')
                options64.Ltol2 = options.UpdateTol;
            else
                %this seems to be the old default value
                options64.Ltol2 = 10;
            end
            if isfield(options,'DropTol')
                options64.small = options.DropTol;
            end
            if isfield(options,'Utol1')
                options64.Utol1 = options.Utol1;
            end
            if isfield(options,'Utol2')
                options64.Utol2 = options.Utol2;
            end
            if isfield(options,'Uspace')
                options64.Uspace = options.Uspace;
            end
            if isfield(options,'Dense1')
                options64.dens1 = options.Dense1;
            end
            if isfield(options,'Dense2')
                options64.dens2 = options.Dense2;
            end
            if isfield(options,'nzinit')
                options64.nzinit = options.nzinit;
            end
        else
            % generate default options from: options = lusol.luset();
            %   |--------+----------+----------------------------------------------------|
            %   | param  |  default | description                                        |
            %   |--------+----------+----------------------------------------------------|
            %
            %   lusol_mex options
            %   |--------+----------+----------------------------------------------------|
            %   | nzinit |        0 | minimum length for storage arrays                  |
            %   | minit  |        0 | minimum length for 'row-space' arrays              |
            %   | ninit  |        0 | minimum length for 'column-space' arrays           |
            %   | minpad |        5 | minimum padding for row and column space arrays    |
            %   |--------+----------+----------------------------------------------------|
            %
            %   LUSOL integer parameters
            %   |--------+----------+----------------------------------------------------|
            %   | maxcol |        5 | max num cols searched for piv element              |
            %   | pivot  |    'TPP' | pivoting method {'TPP','TRP','TCP','TSP'}          |
            %   | keepLU |        1 | keep the nonzeros, if 0, permutations are computed |
            %   |--------+----------+----------------------------------------------------|
            %
            %   LUSOL real parameters
            %   |--------+----------+----------------------------------------------------|
            %   | Ltol1  |     10.0 | max Lij allowed, default depends on pivot method   |
            %   | Ltol2  |     10.0 | max Lij allowed during updates                     |
            %   | small  |  eps^0.8 | absolute tolerance for treating reals as zero      |
            %   | Utol1  | eps^0.67 | absolute tol for flagging small diags of U         |
            %   | Utol2  | eps^0.67 | rel tol for flagging small diags of U              |
            %   | Uspace |      3.0 |                                                    |
            %   | dens1  |      0.3 |                                                    |
            %   | dens2  |      0.5 |                                                    |
            %   |--------+----------+----------------------------------------------------|

            
            mylu=lusol(1);
            options64 = mylu.luset();
            if 0
                options64.Ltol1 = 10;
                options64.Ltol2 = 10;
                options64.pivot = 'TRP';
                %options64.nzinit = 1e7; %increase if low memory leading to segfault
            end
        end
        options=options64;
        mylu = lusol(A,options);
        
        stats = mylu.stats();
        options.Inform = stats.inform;
        options.Rank   = mylu.rank();%changed  May 27th 2011
        options.Nsing  = stats.nsing;
        options.Growth = stats.growth;
        
        
        
        
        
        
        %matrices
        L = mylu.L0();
        U = mylu.U();
        % row permutation
        p = mylu.p();
        % column permutation
        q = mylu.q();
        
        %    lu.factorize inform codes:
        switch stats.inform
            case 0
                if PrintLevel>0
                    fprintf('%s\n','The LU factors were obtained successfully.')
                end
            case 1
                if 0
                    fprintf('%s\n','U appears to be singular, as judged by lu6chk.');
                end
            case 3
                fprintf('%s\n','Some index pair indc(l), indr(l) lies outside the matrix dimensions 1:m , 1:n.');
            case    4
                fprintf('%s\n','Some index pair indc(l), indr(l) duplicates another such pair.');
            case    7
                fprintf('%s\n','The arrays a, indc, indr were not large enough.');
                fprintf('%s\n','Their length "lena" should be increase to at least');
                fprintf('%s\n','The value "minlen" given in luparm(13).');
            case 8
                fprintf('%s\n','There was some other fatal error.  (Shouldn''t happen!)');
            case 9
                fprintf('%s\n','No diagonal pivot could be found with TSP or TDP.');
                fprintf('%s\n','The matrix must not be sufficiently definite or quasi-definite.');
        end

        
        if stats.inform~=0 && stats.inform~=1
            % solve Ax=b
            b = ones(size(A,1),1);
            x = mylu.solve(b);
            % multiply Ax
            b2 = mylu.mulA(x);
            % check the result
            fprintf('%s\t%g\n','Check norm(b-b2) : ', norm(b-b2))
        end
    otherwise
        [L,U,p,q] = lu(A);
end

%this seems superfluous -Ronan
if m > n
    % U(p(n+1:m),:) = [];               % Remove unwanted rows of U
    % Im = speye(m);
    % L  = [L Im(:,n+1:m)];
end

