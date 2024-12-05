function masse=getmasse(PC)

masse = PC.masse;

if isempty(masse)
    error('masse pas calculee : utiliser calc_masse')
end

