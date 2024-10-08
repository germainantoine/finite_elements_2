function [tilde_AA, tilde_LL] = elimine(AA,LL,Refneu,Coorneu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elimine_stokes :
% Mise en oeuvre de la méthode de pseudo-élimination
%
% SYNOPSIS [tilde_AA, tilde_LL] = elimine_stokes(AA,LL,Refneu,Coorneu)
%          
% INPUT *   AA : matrice de la FV à l'intérieur du maillage
%           LL : matrice du second membre sur l'intérieur
%           Refneu : tableau de référence des noeuds
%
% OUTPUT -  tilde_AA : matrice de la FV sur le maillage
%           tilde_LL : matrice du second membre
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=size(Refneu);
%Ns=size(AA(:,1))-2*N;
tilde_AA = AA;
tilde_LL = LL;
for i=1:N
    if Refneu(i)==1 %on sélectionne les noeuds du bord
        tilde_LL(i)= 0;%g(Coorneu(i,1),Coorneu(i,2)); %on applique la condition de Dirichlet sur u
        for j=1:N
            if i==j
                tilde_AA(i,j)=1; %les éléments diagonaux valent 1 sur u
            else
                tilde_AA(i,j)=0; % on élimine les éléments non diagonaux sur u
                tilde_AA(j,i)=0;
            end
        end
    end
end

end
