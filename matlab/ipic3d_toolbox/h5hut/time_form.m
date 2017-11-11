function tf = time_form(t)
global randi1 randi2
T=100;
tf=0;
for k=1:30
    tf= tf +(2*randi1(k)-1)* sin(k.*t/T)+ (2*randi2(k)-1)*cos(k.*t/T);
end
tf./4;