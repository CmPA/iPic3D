%%% all this mess in case I need to read particle info I have extracted
%%% with PythonPars.py into the usual MATLAB format
%%% NB: this works as long as everything is periodic (number of particles stays the same)

clear

results_dir_PY='/Users/innocent/Documents/ExpandingBox/iPic3D/NonEB_Lx10_resonantPars/'
results_dir_PY='/Users/innocent/Documents/ExpandingBox/iPic3D/Resonance_EB_SavingOften_ByBz_EyEz_n1_Lx10_NewINIT_deltaB_002_Longer/'
results_dir_PY='/Users/innocent/Documents/ExpandingBox/iPic3D/NonEB_Lx10_resonantPars_scr_R230_dt4/'
%results_dir_PY='/Users/innocent/Documents/Aurora/Aurora_Max_Lx10_2D/'

ww= [results_dir_PY 'partFromPy/part_proc0.mat']

%%% so all the other saves to workspace are going to be appends
save ([results_dir_PY 'workspace_PyToMat.mat'], 'results_dir_PY' )

load(ww, 'cycles', 'NS','NumHDF')
TOT_HDF=NumHDF;
x_time= cycles';

save ([results_dir_PY 'workspace_PyToMat.mat'], 'x_time' , '-append')
Nsp=NS;

var='qxyzuvw';
str1=[];

% outer cycle should be on particle species
for sp=0:Nsp 

    nameList= [];
    for v=1:size(var,2) 
    
        nameList= [nameList ;  var(v) num2str(sp)];
    end


    nameList

for nl=1:max(size(nameList))
    val= nameList(nl,:)
    str1=[];
    for cyc=1:max(size(cycles))
        for h =0:TOT_HDF-1
            str1= [str1 ' ' val '_pr' num2str(h) '{' num2str(cyc) '}'];
        end
        str1= [str1 ';'];
    end


    
    for h =0:TOT_HDF-1
        ww=  [results_dir_PY 'partFromPy/part_proc' num2str(h) '.mat'];
        load (ww, [val '*']);
    end 
    display(str1);
    str2= [ val '=[' str1 '];'];
    eval(str2)
    
    display(str2) 
    
    %%% traspose, to have the same dimensions as matlab output
    
    display('before eval') 
    eval([val '=' val ''';'])
     
    eval(['clear ' val '_pr* '])
    
    save ([results_dir_PY 'workspace_PyToMat.mat'], val, '-append' ) % this should be faster
    
    eval(['clear ' val ])
    display (['finished ' val])
end
 
end
% end cycle on particle species


%%% if NS>2, combine electron and ion species
if (Nsp>2)
    load([results_dir_PY 'workspace_PyToMat.mat'], 'x0', 'x2' )
    x0=vertcat(x0, x2);
    save([results_dir_PY 'workspace_PyToMat.mat'], 'x0' , '-append' ); clear x0 x2
    
    load([results_dir_PY 'workspace_PyToMat.mat'], 'y0', 'y2' )
    y0=vertcat(y0, y2);
    save([results_dir_PY 'workspace_PyToMat.mat'], 'y0' , '-append'); clear y0 y2
   
    load([results_dir_PY 'workspace_PyToMat.mat'], 'z0', 'z2' )
    z0=vertcat(z0, z2);
    save([results_dir_PY 'workspace_PyToMat.mat'], 'z0' , '-append'); clear z0 z2
    
    load([results_dir_PY 'workspace_PyToMat.mat'], 'u0', 'u2' )
    u0=vertcat(u0, u2);
    save([results_dir_PY 'workspace_PyToMat.mat'], 'u0', '-append' ); clear u0 u2
    
    load([results_dir_PY 'workspace_PyToMat.mat'], 'v0', 'v2' )
    v0=vertcat(v0, v2);
    save([results_dir_PY 'workspace_PyToMat.mat'], 'v0', '-append' ); clear v0 v2
    
    load([results_dir_PY 'workspace_PyToMat.mat'], 'w0', 'w2' )
    w0=vertcat(w0, w2);
    save([results_dir_PY 'workspace_PyToMat.mat'], 'w0', '-append' ); clear w0 w2
    
    load([results_dir_PY 'workspace_PyToMat.mat'], 'q0', 'q2' )
    q0=vertcat(q0, q2);
    save([results_dir_PY 'workspace_PyToMat.mat'], 'q0', '-append' ); clear q0 q2
end

if (Nsp>2)
    load([results_dir_PY 'workspace_PyToMat.mat'], 'x1', 'x3' )
    x1=vertcat(x1, x3);
    save([results_dir_PY 'workspace_PyToMat.mat'], 'x1' , '-append'); clear x1 x3
    
    load([results_dir_PY 'workspace_PyToMat.mat'], 'y1', 'y3' )
    y1=vertcat(y1, y3);
    save([results_dir_PY 'workspace_PyToMat.mat'], 'y1', '-append' ); clear y1 y3
    
    load([results_dir_PY 'workspace_PyToMat.mat'], 'z1', 'z3' )
    z1=vertcat(z1, z3);
    save([results_dir_PY 'workspace_PyToMat.mat'], 'z1', '-append' ); clear z1 z3
    
    load([results_dir_PY 'workspace_PyToMat.mat'], 'u1', 'u3' )
    u1=vertcat(u1, u3);
    save([results_dir_PY 'workspace_PyToMat.mat'], 'u1', '-append' ); clear u1 u3
    
    load([results_dir_PY 'workspace_PyToMat.mat'], 'v1', 'v3' )
    v1=vertcat(v1, v3);
    save([results_dir_PY 'workspace_PyToMat.mat'], 'v1', '-append' ); clear v1 v3
    
    load([results_dir_PY 'workspace_PyToMat.mat'], 'w1', 'w3' )
    w1=vertcat(w1, w3);
    save([results_dir_PY 'workspace_PyToMat.mat'], 'w1', '-append' ); clear w1 w3
    
    load([results_dir_PY 'workspace_PyToMat.mat'], 'q1', 'q3' )
    q1=vertcat(q1, q3);
    save([results_dir_PY 'workspace_PyToMat.mat'], 'q1', '-append' ); clear q1 q3
end