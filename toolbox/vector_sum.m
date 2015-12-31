function v = vector_sum(v1, v2, n)

v = zeros(n, 1);

if length(v1) < n
    for i = 1 : length(v1)
        v(i) = v(i) + v1(i);
    end
elseif length(v1) > n
    v = v + v1(1:n);
else
    v = v + v1;
end

if length(v2) < n
    for i = 1 : length(v2)
        v(i) = v(i) + v2(i);
    end
elseif length(v2) > n
    v = v + v2(1:n);
else
    v = v + v2;
end


