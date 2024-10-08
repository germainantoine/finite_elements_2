function [tilde_AA, tilde_LL] = elimine_stokes(AA,LL,Refneu,Coorneu,nbpt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elimine_stokes :
% Mise en oeuvre de la méthode de pseudo-élimination
%
% SYNOPSIS [tilde_AA, tilde_LL] = elimine_stokes(AA,LL,Refneu,Coorneu)
%          
% INPUT *   AA : matrice de la FV à l'intérieur du maillage
%           LL : matrice du second membre sur l'intérieur
%           Refneu : tableau de référence des noeuds
%           Coorneu : coordonnées des noeuds du maillage
%
% OUTPUT -  tilde_AA : matrice de la FV sur le maillage
%           tilde_LL : matrice du second membre
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=nbpt;

tilde_AA = AA;
tilde_LL = LL;
for i=1:N
    if Refneu(i)==1
        tilde_LL(i)=g1(Coorneu(i,1),Coorneu(i,2));
        tilde_LL(i+N)=g2(Coorneu(i,1),Coorneu(i,2));
        tilde_AA(i,:)=0;
        tilde_AA(i+N,:)=0;
        tilde_AA(i,i)=1;
        tilde_AA(N+i,N+i)=1;
    end
        if Refneu(i)==2
            tilde_LL(i)=g1(Coorneu(i,1),Coorneu(i,2));
            tilde_LL(i+N)=g2(Coorneu(i,1),Coorneu(i,2));
            tilde_AA(i,:)=0;
            tilde_AA(i+N,:)=0;
            tilde_AA(i,i)=1;
            tilde_AA(N+i,N+i)=1;
        end
    end

end