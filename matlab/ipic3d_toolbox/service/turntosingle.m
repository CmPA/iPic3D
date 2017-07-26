variable=whos;
[nvar dummy]=size(variable);
for ivar=1:nvar
if eval(['isreal(' variable(ivar).name ')'])
eval([variable(ivar).name '=' 'single(' variable(ivar).name ');']);
%disp([num2str(ivar), variable(ivar).name]);
end
end

