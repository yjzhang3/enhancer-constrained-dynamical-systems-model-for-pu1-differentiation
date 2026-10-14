function p = parentdir(path, N)
    for k = 1:N+1
        path = fileparts(path);
    end
    p = path;
end