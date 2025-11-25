function [z] = rectify(x)
z = x .*(x>0); 
end