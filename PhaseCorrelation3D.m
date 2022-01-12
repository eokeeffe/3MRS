function [ m, n] = PhaseCorrelation3D( ref, sen, matchrad)

fftRef = fftn(ref);
fftSen = fftn(sen);
fftRef_conj = conj(fftRef);
corr = fftSen.*fftRef_conj;    
corr = real(fftshift(ifftn( corr )));

max_corr = max( corr( : ) );
max_index = find( corr == max_corr );
[ m, n, ~ ] = ind2sub( size( corr ), max_index );

m = m - 1 - matchrad;
n = n - 1 - matchrad;