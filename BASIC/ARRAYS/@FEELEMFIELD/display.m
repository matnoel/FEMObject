function display(u)
% function display(u)

disp(class(u))
for k=1:length(u.value)
    disp(['Element group #' num2str(k)])
    disp(u.value{k})
end