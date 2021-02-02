name_elmt = 'mesh2d_ex/mesh.elements';
name_nodes = 'mesh2d_ex/mesh.nodes';



%import elements and nodes files
elmt = importdata(name_elmt);
nodes = importdata(name_nodes);



%check each element
cmp_degen = 0;
cmp_align = 0;
for k=1: length(elmt)
    if elmt(k,3) ~= 706
        disp('Code only developped for 706 elements');
        disp(['Here element number', num2str(k), 'is of type' num2str(elmt(k,3))]);
        return;
    else
        % nodes of the elements
        n1 = elmt(k,4);
        n2 = elmt(k,5);
        n3 = elmt(k,6);
        n4 = elmt(k,7);
        n5 = elmt(k,8);
        n6 = elmt(k,9);
        %coordinate of the nodes
        xn1 = nodes(n1, 3);
        yn1 = nodes(n1, 4);
        zn1 = nodes(n1, 5);
        xn2 = nodes(n2, 3);
        yn2 = nodes(n2, 4);
        zn2 = nodes(n2, 5);
        xn3 = nodes(n3, 3);
        yn3 = nodes(n3, 4);
        zn3 = nodes(n3, 5);
        xn4 = nodes(n4, 3);
        yn4 = nodes(n4, 4);
        zn4 = nodes(n4, 5);
        xn5 = nodes(n5, 3);
        yn5 = nodes(n5, 4);
        zn5 = nodes(n5, 5);
        xn6 = nodes(n6, 3);
        yn6 = nodes(n6, 4);
        zn6 = nodes(n6, 5);
        %check vertical alignment
        if xn1~=xn4 || yn1~=yn4 || xn2~=xn5 || yn2~=yn5 || xn3~=xn6 || yn3~=yn6
            disp('_________________________________________');
            disp(['nodes not vertically aligned for elents: ', num2str(k)]);
            cmp_align = cmp_align + 1;
        end
        if zn4-zn1<=0 || zn5-zn2<=0 || zn6-zn3<=0
            disp('_________________________________________');
            disp(['Element degenerated: ', num2str(k)]);
            disp(['x1, y1, z1, z4: ' num2str(xn1) ', ' num2str(yn1) ', ' num2str(zn1) ', ' num2str(zn4)]); 
            disp(['x2, y2, z2, z5: ' num2str(xn2) ', ' num2str(yn2) ', ' num2str(zn2) ', ' num2str(zn5)]); 
            disp(['x3, y3, z3, z6: ' num2str(xn3) ', ' num2str(yn3) ', ' num2str(zn3) ', ' num2str(zn6)]); 
            cmp_degen = cmp_degen + 1;
            
        end
            
        
       
        
        
        
        
    end
end

disp('_____________________________');
disp([num2str(cmp_degen) 'degenerated elements over ' num2str(length(elmt)) ' elements']);

disp('_____________________________');
disp([num2str(cmp_align) ' not vertical aligned elements over ' num2str(length(elmt)) ' elements']);























