function setglob(name,val)
eval(['global ' name]);
eval([name '.field1=' num2str(val)])
