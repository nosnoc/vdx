warning("Installing vdx not through nosnoc. If you call `install_nosnoc` later then this may cause conflicts");

[vdx_path,~,~] = fileparts(mfilename("fullpath"));
sep = pathsep;

root_already_installed = contains([sep, path, sep], [sep, vdx_path, sep], 'IgnoreCase', ispc)

if ~root_already_installed
    addpath(nosnoc_path);
end
savepath
