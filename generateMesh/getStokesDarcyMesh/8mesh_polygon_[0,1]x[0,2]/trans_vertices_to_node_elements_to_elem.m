for n = 1:6
    meshname = 'poly';
    Smeshname = ['Smesh_',meshname,'_[0,1]x[1,2]_',num2str(2^(n))];
    Dmeshname = ['Dmesh_',meshname,'_[0,1]x[0,1]_',num2str(2^(n))];
    
    mesh_name = ['mesh_Darcy_',num2str((2^(n))^2)];
    load(mesh_name);
    node = vertices';
    elem = elements';
    save(Dmeshname,'node','elem'); clear node elem
    
    mesh_name = ['mesh_Stokes_',num2str((2^(n))^2)];
    load(mesh_name);
    node = vertices';
    elem = elements';
    save(Smeshname,'node','elem');
end