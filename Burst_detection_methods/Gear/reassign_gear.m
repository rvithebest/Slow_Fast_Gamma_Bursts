function params_out = reassign_gear(params_in,R,sRate)
L = length(R);
R = R';
t_p1 =params_in(1);
f_p1 = params_in(2);
s_p1 = params_in(3);
t_delta=0.1/L;
f_delta = 0.1/sRate;
t_p11 = t_p1+t_delta;
t_p1_1 = t_p1-t_delta;
f_p11 = f_p1+f_delta;
f_p1_1 = f_p1-f_delta;

g_p1 = mageAtom_real([t_p1 f_p1 s_p1],sRate,L,0)+1i*mageAtom_imag([t_p1 f_p1 s_p1],sRate,L,0);
g_p11 = mageAtom_real([t_p11 f_p11 s_p1],sRate,L,0)+1i*mageAtom_imag([t_p11 f_p11 s_p1],sRate,L,0);
g_p1_1 = mageAtom_real([t_p1_1 f_p1_1 s_p1],sRate,L,0)+1i*mageAtom_imag([t_p1_1 f_p1_1 s_p1],sRate,L,0);
g_p12 = mageAtom_real([t_p11 f_p1_1 s_p1],sRate,L,0)+1i*mageAtom_imag([t_p11 f_p1_1 s_p1],sRate,L,0);

F = [((f_p11^2-f_p1^2)-exp(-2*s_p1)*(t_p11^2-t_p1^2)) (t_p1-t_p11) (f_p1-f_p11);((f_p1_1^2-f_p1^2)-exp(-2*s_p1)*(t_p1_1^2-t_p1^2)) (t_p1-t_p1_1) (f_p1-f_p1_1);((f_p1_1^2-f_p1^2)-exp(-2*s_p1)*(t_p11^2-t_p1^2)) (t_p1-t_p11) (f_p1-f_p1_1)];
A = [(log(abs(R*g_p1')/abs(R*g_p11'))-(pi*exp(-s_p1)*(t_p11^2-t_p1^2)));(log(abs(R*g_p1')/abs(R*g_p1_1'))-(pi*exp(-s_p1)*(t_p1_1^2-t_p1^2)));(log(abs(R*g_p1')/abs(R*g_p12'))-(pi*exp(-s_p1)*(t_p11^2-t_p1^2)))];
B = F\A;
b = B(1);
B_t = B(2);
B_f = B(3);

f_t = B_f/(2*b);
t_t = B_t/(2*(exp(-2*s_p1))*(pi*exp(s_p1)-b));

s_t = s_p1-log(pi/b*exp(s_p1)-1);
if isreal(s_t)
    params_out = [t_t f_t s_t];
else
    params_out = params_in;
end



end