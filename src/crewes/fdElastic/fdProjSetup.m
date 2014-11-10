function modelDir = fdProjSetup (project_dir)
if (nargin < 1)
    project_dir = 0;
end
if (project_dir == 0)
    modelDir = uigetdir([],'Select mFD2D project directory');
    if (modelDir == 0)
       warning('mFD2D:cancelled','User cancelled project directory selection');
       return;
    else
        project_dir = modelDir;
    end
else
    modelDir = project_dir;
end
disp(project_dir)
