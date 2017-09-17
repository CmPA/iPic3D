function  divA = compute_div(x,y,z,Ax,Ay,Az)
divA = divergence(x,y,z,permute(Ax,[2 1 3]), permute(Ay, [2 1 3]), permute(Az, [2,1,3]));
divA = permute(divA, [2 1 3]);
end