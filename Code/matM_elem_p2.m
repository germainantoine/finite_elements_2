function [Mel] = matM_elem_p2(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de la matrice de masse elementaire en P2 Lagrange
%
% SYNOPSIS [Mel] = matM_elem_p2(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - Mel matrice de masse elementaire (matrice 6x6)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% calcul de la matrice de masse
% -------------------------------
% Initialisation
Mel = zeros(6); % A COMPLETER

% Points et poids de quadrature
S_hat = [0.0915762135098, 0.0915762135098;
         0.8168475729805, 0.0915762135098;
         0.0915762135098, 0.8168475729805;
         0.1081030181681, 0.4459484909160;
         0.4459484909160, 0.1081030181681;
         0.4459484909160, 0.4459484909160];
poids = [0.05497587183, 0.05497587183, 0.05497587183, 0.1116907948, 0.1116907948, 0.1116907948];

% A COMPLETER
determinantBl = abs((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));
for i=1:6
    for j=1:6
        for q=1:6
            wq=[(-S_hat(q,1)-S_hat(q,2)+1)*(-2*S_hat(q,1)-2*S_hat(q,2)+1),
                S_hat(q,1)*(2*S_hat(q,1)-1),
                S_hat(q,2)*(2*S_hat(q,2)-1),
                4*(1-S_hat(q,1)-S_hat(q,2))*S_hat(q,1),
                4*S_hat(q,1)*S_hat(q,2),
                4*S_hat(q,2)*(1-S_hat(q,1)-S_hat(q,2))];
            Mel(i,j) = Mel(i,j) + poids(q)*wq(i)*wq(j)*determinantBl;
        end
    end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2023
