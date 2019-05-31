%
%   This is a simulation of the ChR2 masking experiment in 
%   "Concentration invariant odor coding "
%   by Christopher D. Wilson, Gabriela O. Serrano, Alexei A. Koulakov, Dmitry Rinberg
%   The simulation includes olfactory bulb-to-piriform cortex network (PC) and 
%   PC-to-PC network in the firing rate model.
%   Questions: akula@cshl.edu (Alex Koulakov)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all


%load masks
%load mc_ord

%Nan = 4;            % number of animals
Ntr = 100;




T_chr2 = [0.05:0.05:1];
N_chr2 = length(T_chr2);

Smells = [zeros(Ntr/2,N_chr2);ones(Ntr/2,N_chr2)];           % 0=A, 1=B
Resp = Smells*0;
OA = Resp;
OB = Resp;

%[ f, p, F, P, B ] = run1sim1( 'B', [1:Ng]', [8, 6, 1, 5, 3, 7, 4, 2, 9:Ng]', Ng, Np, Nt, dt, .05, thU, thL, 0, 2, chr2, t_sm, PBExFwd, WPP*0.6, WePB*0.025, cMult, 1000, 1 );

%
%   Parameters of simulation
%

%t_chr2 = 0.7;                       % time of chr2 onset

thL = -150;         % lower threshold 
thU = 0.2;            % upper threshold % 0.1-

Ng = 300; 
Np = 1000; %

NePP = 3;
NiPP = 30;
NePB = 40;
NeSelf = 3;
NiSelf = 500;

PBExFwd = 6*0.025;
PPExSelf = 5*0.6;
ANoise = 0.1;
%chr2Amp=0.075;                        % amplitude of chr2 stim
chr2Amp=0.18;                        % amplitude of chr2 stim
cMult = 1;
cAdd = 100.5;                       % shift added between two sequences that determines mixture conc
                                    % set to 2.5 for really complex misture
TransientP = 0.9;                   % probability of a transient to be emmitted. It is concentration dependent. 
                                    % choose in the range 0.8-0.9


PPIn = -5*0.6; %-0.65;

tau = 0.05;
dt = 0.002;
Nt = round(1.2/dt);

t_sm=0.4;                           % time when odorant information is delivered







if 0
    %s2 = RandStream.create('mt19937ar','Seed',2466);       
    %prevstream = RandStream.setGlobalStream(s2);

    MMM = randperm(Ng)';
    MMM = randperm(Ng)';

    for i=1:Ng
        dd = randn; 
    end
end




for n_chr2 = [1: N_chr2]
    n_chr2
    tic

    
     t_chr2 = T_chr2(n_chr2);

    
     rrr = zeros(size(Resp(:,n_chr2)));
     
     %parfor (ntr = 1:Ntr, 4)
     
     for ntr = 1:Ntr
        
         if  mod(ntr, 10)==1                % this will reset all the variables 
                                            % initiaing a new animal
                                            
                                            %ntr
             
             % Setup new weights to mimic a new animal
             
             MC_ordA = randperm(Ng)';
             MC_ordB = randperm(Ng)';
             
             WePB = random_matrix1(Np, Ng, NePB);           % weights bulb -> piriform
             WePP = random_matrix1(Np,Np,NeSelf)*PPExSelf; 
             WiPP = random_matrix1(Np,Np,NiSelf)*PPIn; 
             WPP = WePP+WiPP;
             WPP = WPP*0.25;
             WePB = random_matrix1(Np, Ng, NePB)*0.25*PBExFwd;

             chr2 = chr2Amp*(rand(Ng,1) > 0.3);

            [ maskA, p, F, P, B ] = run1sim2( 'A', MC_ordA, MC_ordB, Ng, Np, Nt, dt, tau, thU, thL, 0, 1, 2, chr2, t_sm, WPP, WePB, cMult, cAdd, 0 );
            [ maskB, p, F, P, B ] = run1sim2( 'B', MC_ordA, MC_ordB, Ng, Np, Nt, dt, tau, thU, thL, 0, 1, 2, chr2, t_sm, WPP, WePB, cMult, cAdd, 0 );

         end
         
        if Smells(ntr,n_chr2)
            Smell='B';
        else
            Smell='A';
        end


        [ f, p, F, P, B ] = run1sim2( Smell, MC_ordA, MC_ordB, Ng, Np, Nt, dt, tau, thU, thL, ANoise, TransientP, t_chr2, chr2, t_sm, WPP, WePB, cMult, cAdd, 0 );


        ff = full(f);
        oA = sum(maskA.*ff);
        oB = sum(maskB.*ff);

        OA(ntr, n_chr2) = oA;
        OB(ntr, n_chr2) = oB;

        smell_presented = Smell;

        if oA>oB
            smell_detected = 'A';
        else
            smell_detected = 'B';
        end

        rrr(ntr)=(smell_detected=='B');

    end, 



Resp(:, n_chr2) = rrr ;



toc
Correct = (Smells == Resp);


figure(2)
subplot(2,2,2)
imagesc(OA), colorbar

subplot(2,2,4)
imagesc(OB), colorbar


if 0
subplot(2,2,1)
D1 = zeros(25,40,3);
D1(:,:,1) = reshape(f,25,40);
D1(:,:,2) = reshape(maskA,25,40);
image(D1), axis equal, axis off
colorbar

subplot(2,2,3)
D1(:,:,2) = reshape(maskB,25,40);
image(D1), axis equal, axis off
colorbar
end

drawnow






figure(1)
plot(T_chr2, mean(Correct), '.-'), drawnow

end






