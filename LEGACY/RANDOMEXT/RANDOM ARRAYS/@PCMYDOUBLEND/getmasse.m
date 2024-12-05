function m = getmasse(u)

if isradial(u)
    m = getximasse(u.L); 
else
    m = getmasse(u.L);    
end