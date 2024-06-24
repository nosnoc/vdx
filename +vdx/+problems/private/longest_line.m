function n = longest_line(str)
    idx = [regexp(str, '\\n'), length(str)];
    if length(idx) == 1
        n = idx;
    else
        n = max(diff(idx))-1;
    end
end
