function y = callback_localfft(x,dir,options)


if isfield(options, 'w')
    w = options.w;
else
    w = 9;
end
w = 9;
if isfield(options, 'q')
    q = options.q;
else
    q = round((w-1)/2);
end
if isfield(options, 'n')
    n = options.n;
else
    error('You must specify options.n');
end

if dir==1
    y = perform_windowed_fourier_transform(x,w,q,n, options);
else
    y  = perform_windowed_fourier_transform(x,w,q,n, options);
end