% %---------------------------------------------------------------------------------------------------------------------
% octave m file
%
% Author:		Jan Stegenga
% Company:		INCAS3
% Date:			3/10/2014
% Copyright:	INCAS3 2016
% License:		GNU GPLv3.0
% Requirements:	32bit octave (i used windows 7).
% 
% Description:
%	Illustrates the cyclic reduction method of solving b = A*x. 
%	where A is a diagonal matrix with 3 diagonals, as often found in (discretized) spatio-temporal dynamic systems. 
%	Here, size( A ) = [N, N], where N is a power of 2.
% %--------------------------------------------------------------------------------------------------------------------


% %--------------------------------------------------------------------------------------------------------------------
% shifts vector d by k to left or right and inserts k element of value <num>. d is expected to be a row vector.
% %--------------------------------------------------------------------------------------------------------------------
function out = shiftvector( d, k, num )
	if k<0
		out = [num*ones(1,abs(k)), d](1:end+k);
	elseif k>0
		out = [d, num*ones(1,k)](k+1:end);
	else
		out = d;
	end
endfunction

% %--------------------------------------------------------------------------------------------------------------------
% generate the constants used for cyclic reduction. a, b, c, are the above, on and below diagonal elements of matrix A. 
% size(A) = [N,N] elements, where N is a power of 2.
% %--------------------------------------------------------------------------------------------------------------------
function [F,G,H] = generate_cyclic_reduction_vectors( a, b, c, N )
	steps = log2(N)+1;
	ank = zeros( N, steps );
	bnk = zeros( N, steps );
	cnk = zeros( N, steps );
	ank( :, 1 ) = [0;a(:)];
	bnk( :, 1 ) = b(:);
	cnk( :, 1 ) = [c(:);0];
	for k=2:steps
		ank( :, k ) =   - ank( :, k-1 ).*shiftvector( ank( :, k-1 )', -2^(k-2), 0 )'./shiftvector( bnk( :, k-1 )', -2^(k-2), 1 )';
		bnk( :, k ) = 	  bnk( :, k-1 ) ...
						- ank( :, k-1 ).*shiftvector( cnk( :, k-1 )', -2^(k-2), 0 )'./shiftvector( bnk( :, k-1 )', -2^(k-2), 1 )' ...
						- cnk( :, k-1 ).*shiftvector( ank( :, k-1 )', +2^(k-2), 0 )'./shiftvector( bnk( :, k-1 )', +2^(k-2), 1 )';
		cnk( :, k ) =   - cnk( :, k-1 ).*shiftvector( cnk( :, k-1 )', +2^(k-2), 0 )'./shiftvector( bnk( :, k-1 )', +2^(k-2), 1 )';
	end

	for k=2:steps
		F( :, k ) = ank( :, k-1 )./shiftvector( bnk( :, k-1 )', -2^(k-2), 1 )'; 
		G( :, k ) = cnk( :, k-1 )./shiftvector( bnk( :, k-1 )', +2^(k-2), 1 )';
	end
	H = 1./bnk( :, steps );
endfunction

% %--------------------------------------------------------------------------------------------------------------------
% perform the cyclic reduction 
% %--------------------------------------------------------------------------------------------------------------------
function x_by_cr = cyclic_reduction( F, G, H, b )
	x_by_cr = b;
	for k = 2:( log2(length(b)) + 1 )
		x_by_cr = x_by_cr -F(:,k) .* shiftvector( x_by_cr', -2^(k-2), 0 )' -G(:,k) .* shiftvector( x_by_cr', +2^(k-2), 0 )';	
	end
	x_by_cr = x_by_cr .* H;
endfunction


%-----------
%usage:
% find the solution to x for the matrix equation: b = A*x,
% where A is of the form: 
% b		c		0		0		...		0		0		0
% a		b		c		0		...		0		0		0
% 0		a		b		c		...		0		0		0
% ...
% 0		0		0		0				a		b		c
% 0		0		0		0				0		b		c
% 0		0		0		0				0		0		b
%
%A general solution is to use gaussian elimination, but that offers limited possibilities for parallelization of calculations.
%-----------

% 16x16 example
N = 16;

% matrix A is diagonal [-1, 0, +1] 
aA = ones(N-1, 1);											% N above-diagonal elements
bA = [ - 1; - 2.5*ones( N-2,1 );- 4  ];					% N+1 diagonal elements
cA = ones(N-1, 1);											% N below-diagonal elements
A = diag( bA, 0 ) + diag( aA, 1 ) + diag( cA, -1);							

% use random numbers for b
b = randn( N, 1 );

%one solution is to use the inverse of A, or the '\' operator in many languages which uses efficient BLAS operations to find the solution
invA = inv(A);
x_by_inverse = invA*b; 

% construct the matrices used for cyclic reduction
[F, G, H] 	= generate_cyclic_reduction_vectors( aA, bA, cA, N );

% solve the matrix equation
x_by_cr 	= cyclic_reduction( F, G, H, b ); 

% compare results
[x_by_inverse, x_by_cr] 




