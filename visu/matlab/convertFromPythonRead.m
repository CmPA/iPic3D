%%% all this mess in case I need to read particle info I have extracted
%%% with PythonPars.py into the usual MATLAB format
%%% NB: this works as long as everything is periodic (number of particles stays the same)

clear

results_dir_PY='/Users/innocent/Documents/ExpandingBox/iPic3D/NonEB_Lx10_resonantPars/'
results_dir_PY='/Users/innocent/Documents/ExpandingBox/iPic3D/Resonance_EB_SavingOften_ByBz_EyEz_n1_Lx10_NewINIT_deltaB_002_Longer/'
results_dir_PY='/Users/innocent/Documents/ExpandingBox/iPic3D/NonEB_Lx10_resonantPars_scr_R230_dt4/'
%results_dir_PY='/Users/innocent/Documents/Aurora/Aurora_Max_Lx10_2D/'

ww= [results_dir_PY 'partFromPy/part_proc0.mat']

load(ww)
TOT_HDF=NumHDF;
x_time= cycles;
Nsp=NS;

var='qxyzuvw';
str1=[];

nameList= [];
for v=1:size(var,2) 
    for sp=0:2 :Nsp  %only electrons
        nameList= [nameList ;  var(v) num2str(sp)];
    end
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
    
    display (['finished ' val])
end

x_time= x_time';


%%% if NS>2, combine electron and ion species
if (exist('x2'))
    x0=vertcat(x0, x2);
    y0=vertcat(y0, y2);
    z0=vertcat(z0, z2);
    
    u0=vertcat(u0, u2);
    v0=vertcat(v0, v2);
    w0=vertcat(w0, w2);
    
    q0=vertcat(q0, q2);
end

if (exist('x3'))
    x1=vertcat(x1, x3);
    y1=vertcat(y1, y3);
    z1=vertcat(z1, z3);
    
    u1=vertcat(u1, u3);
    v1=vertcat(v1, v3);
    w1=vertcat(w1, w3);
    
    q1=vertcat(q1, q3);
end

clear *2 *3 *pr* va* cyc h nameList N* sp nl  TOT* ww str*
display('saving workspace')
%save ([results_dir_PY 'workspace_PyToMat.mat'], '-v7.3')
save ([results_dir_PY 'workspace_PyToMat.mat'], '-v6') % this should be faster
