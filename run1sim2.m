function [ f, p, F, P, B ] = run1sim1( Smell, MC_ordA, MC_ordB, Ng, Np, Nt, dt, tau, thU, thL, ANoise, TransientP, t_chr2, chr2, t_sm, WPP, WePB, cMult, cAdd, display )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% [ f, p, F, P, B ] = run1sim1( 'B', [1:Ng]', [8, 6, 1, 5, 3, 7, 4, 2, 9:Ng]', Ng, Np, Nt, dt, .05, thU, thL, 0, 2, chr2, t_sm, PBExFwd, WPP*0.6, WePB*0.025, cMult, 1000, 1 );


if display
        figure
end



if 0


    if Smell == 'A'
        MC_ord = MC_ordA;
    else
        MC_ord = MC_ordB;
    end

else
    
    MC_ord = [MC_ordA; MC_ordB];
    ttt = [1:(length(MC_ordA))]';
    
    if Smell == 'A'
        TTT = [ttt; ttt*cMult+cAdd];
    else
        TTT = [ttt*cMult+cAdd; ttt];
    end
    
    [sTTT, ind] = sort(TTT, 'ascend');
    
    MC_ord = MC_ord(ind);
    
    % Remove duplicates
    
    mmm = zeros(size(MC_ordA));
    
    n_inc = 1;
    
    if display
        MC_ord(1:20)'
    end
    
    for i=1:length(MC_ord)
        if MC_ord(i) == 0
            continue
        end
        
        ind = find(MC_ord == MC_ord(i));
        mmm(n_inc) = MC_ord(i);
        n_inc = n_inc + 1;
        MC_ord(ind)=0;
        
    end
    
    MC_ord = mmm;
    
    
end

B = zeros(Ng,Nt);

%dn_bulb = Nt/Ng*1.2;
%dn_bulb = Nt/Ng*2.1;       % this was used in simulation on the night of
%5/22-23/2016


dn_bulb = round(0.2/10/dt);      % 10 primary glomeruli in 0.2 sec


%dn_pulse = round(0.16/dt);
dn_pulse = round(0.5/dt);
n0_pulse = round(t_sm/dt);

if display
    Mitral_cell_order = MC_ord(1:20)'
end

for i=1:Ng
   n_glom = MC_ord(i);
   
   dd=round(randn*0/dt);                      % jitter 
   t1 = round(dn_bulb*i+n0_pulse+dd);
   t2 = round(dn_bulb*i+n0_pulse+dn_pulse+dd);
   t2 = min(t2,Nt);
   
   if t1<Nt
       if rand<TransientP
            B(n_glom,t1:t2)=1;
       end
   end
    
end

NN = ANoise*randn(size(B));
B = B+NN;



I0 = 0*ones(Np,1);


% begin sim...............................................................

act = 0;

p = zeros(Np,1);
f = p;
f_prev = 0;
b = zeros(Ng,1);
p_prev = p;

P = zeros(Np,Nt);
F = P;



step=1;
t=0;


% add chr stim

%for i=round(t_chr2/dt):size(B,2)

for i=round(t_chr2/dt):min(round((t_chr2+0.1)/dt),size(B,2))

    B(:,i) = B(:,i)+chr2;

end

B = sparse(B);
f = sparse(f);
b0 = sparse(B(:,1)*0);

dN = round(0.2/dt);

for i=(-dN):Nt

    if i>0
        b = B(:,i);
    else
        b = ANoise*randn(size(B(:,1)));
    end
    p_prev = p;
    f_prev = f;
        
    t=t+dt;
    step = step+1;
        
    I = I0 + WPP*f + WePB*b;
        
    p = p_prev + dt/tau * (-p_prev + I);
        
    ind = find(p>=thU);
    f(ind)=1;
        
    ind = find(p<=thL);
    f(ind)=0;

    ind = find((thL < p).*(p < thU));
    f(ind) = f_prev(ind);
    
    if display
        
        
        subplot(1,2,1)
    
        image(64*reshape(b,15,20)), axis equal, axis off, colormap winter

        subplot(1,2,2)
    
        image(64*reshape(f,25,40)), axis equal, axis off, colormap winter
        
        title(sprintf('t=%g',i/Nt))
        drawnow
    
        pause(0.02)
    
    end
    
    if i>0
    P(:,i)=p;
    F(:,i)=f;
    end
    
end







end

