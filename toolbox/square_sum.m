function s = square_sum(s1, s2, n)

s = zeros(n, n);

if size(s1, 1) < n
    for i = 1 : size(s1, 1)
        for j = 1 : size(s1, 1)
            s(i, j) = s(i, j) + s1(i, j);
        end
    end
elseif size(s1, 1) > n
    s = s + s1(1:n, 1:n);
else
    s = s + s1;
end

if size(s2, 1) < n
    for i = 1 : size(s2, 1)
        for j = 1 : size(s2, 1)
            s(i, j) = s(i, j) + s2(i, j);
        end
    end
elseif size(s2, 1) > n
    s = s + s2(1:n, 1:n);
else
    s = s + s2;
end