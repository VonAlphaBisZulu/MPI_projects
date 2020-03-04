A = [1 1; 0 2];
b = [2 3];
c = [1 2];
p = [];
X = 0:0.1:15;
Y = 0:0.1:15;

% Primales Problem x=[i j]
for j=Y
    for i=X
        if all(A*[i j]'<=b')
            p(i==X,j==Y) = c*[i j]';
        else
            p(i==X,j==Y) = NaN;
        end
    end
end
X = X(1:size(p,2));
Y = Y(1:size(p,1));
mesh(X,Y,p);

% Duales Problem y=[i j]
p = [];
X = 0:0.1:3;
Y = 0:0.1:3;
for j=Y
    for i=X
        if all(A'*[i j]'>=c')
            p(i==X,j==Y) = b*[i j]';
        else
            p(i==X,j==Y) = NaN;
        end
    end
end
X = X(1:size(p,2));
Y = Y(1:size(p,1));
figure
mesh(X,Y,p);