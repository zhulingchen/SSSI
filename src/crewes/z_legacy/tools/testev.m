function a = testev(name)
eval(['global ' name]);
eval ([name '=5;']);
