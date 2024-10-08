function [Kel] = matK_elem_p2(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de la matrice de rigidité elementaire en P2 Lagrange
%
% SYNOPSIS [Kel] = matK_elem_p2(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - Kel matrice de rigidité elementaire (matrice 6x6)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% calcul de la matrice de rigidité
% -------------------------------
% Initialisation
Kel = zeros(6,6); % A COMPLETER

% Points et poids de quadrature
SCHAPO = [1/6, 1/6;
         2/3, 1/6;
         1/6, 2/3];
poids = 1/6;

% A COMPLETER
Bl=[x2-x1,x3-x1;y2-y1,y3-y1];
det(Bl);
detBl = abs((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));
invTransBl = inv(Bl');

for i=1:6
    for j=1:6
        for q=1:3
            gradwq=[4*SCHAPO(q,1)+4*SCHAPO(q,2)-3 , 4*SCHAPO(q,1)+4*SCHAPO(q,2)-3;
                    4*SCHAPO(q,1)-1, 0;
                    0 , -1+4*SCHAPO(q,2);
                    -8*SCHAPO(q,1)-4*SCHAPO(q,2)+4 , -4*SCHAPO(q,1); 
                    4*SCHAPO(q,2) , 4*SCHAPO(q,1);
                    -4*SCHAPO(q,2) , -4*SCHAPO(q,1)-8*SCHAPO(q,2)+4];
                A=invTransBl;
                B=gradwq(i,:)';
                C=gradwq(j,:)';
                D=A*B;
                E=A*C;
                F=D'*E;
            Kel(i,j) = Kel(i,j) + poids*F*detBl;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2023
