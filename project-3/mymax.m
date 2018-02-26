
function [val] = mymax(a,b,c)
abs_a = abs(a);
abs_b = abs(b);
abs_c = abs(c);


if abs_a>abs_b
    if abs_a>abs_c
        val = a;
    else
        val = c;
    end
else
    if abs_b>abs_c
        val = b;
    else
        val = c;
    end
end

end