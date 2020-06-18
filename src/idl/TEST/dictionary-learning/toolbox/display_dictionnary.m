function H = display_dictionnary(D,X,nb, s)

% H = display_dictionnary(D,X,nb, s);

if nargin<2
    X = [];
end
if nargin<4
    s = 1;
end

if nargin<3
    nb = 10;
end
if length(nb)==1
    nb = [nb nb];
end

w1 = sqrt(size(D,1)/s);
pad = 2;

col = [1 1 1];

if ~isempty(X)
    d = sum( abs(X') );
    [tmp,I] = sort(d);
else
    I = randperm(size(D,2));
end
% nb = min(15, floor(sqrt(K)));
H = repmat( reshape(col,[1 1 3]), [nb*(w1+pad) 1] );
k = 0;
for i=1:nb(1)
    for j=1:nb(2)
        k = k+1;
        if k<length(I)
        v = rescale( D(:,I(k)) );
        v = reshape(v,w1,w1,s);
        selx = (i-1)*(w1+pad)+1:(i-1)*(w1+pad)+w1;
        sely = (j-1)*(w1+pad)+1:(j-1)*(w1+pad)+w1;
        if s==1
            v = repmat(v,[1 1 3]);
        end
        H(selx,sely,:) = v;
        end
    end
end

H(end-pad+1:end,:,:) = [];
H(:,end-pad+1:end,:) = [];

if nargout==0
    clf;
    imagesc(H); axis image; axis off;
end
