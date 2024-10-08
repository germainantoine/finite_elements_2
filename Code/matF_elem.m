function [Fel] = matF_elem(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de la matrice elementaire du bloc rectangulaire (p, dv2/dy)
%
% SYNOPSIS [Fel] = matF_elem(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle
%                      (vecteurs reels 1x2)
%
% OUTPUT - Fel matrice elementaire rectangulaire (matrice 6x3)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% calcul de la matrice elementaire du bloc rectangulaire (p, dv2/dy)
S_hat = [0.0915762135098, 0.0915762135098;
         0.8168475729805, 0.0915762135098;
         0.0915762135098, 0.8168475729805;
         0.1081030181681, 0.4459484909160;
         0.4459484909160, 0.1081030181681;
         0.4459484909160, 0.4459484909160];
poids = [0.05497587183, 0.05497587183, 0.05497587183, 0.1116907948, 0.1116907948, 0.1116907948];

% ------------------------------------------------------------------
Fel = zeros(6,3);
Bl=[x2-x1,x3-x1;y2-y1,y3-y1];
detBl = det(Bl);
invTransBl = inv(Bl');

for i=1:6
    for j=1:3
        for q=1:6
            gradwq=[4*S_hat(q,1)+4*S_hat(q,2)-3 , 4*S_hat(q,1)+4*S_hat(q,2)-3;
                    4*S_hat(q,1)-1, 0;
                    0 , -1+4*S_hat(q,2);
                    -8*S_hat(q,1)-4*S_hat(q,2)+4 , -4*S_hat(q,1); 
                    4*S_hat(q,2) , 4*S_hat(q,1);
                    -4*S_hat(q,2) , -4*S_hat(q,1)-8*S_hat(q,2)+4];
            lambda=[1-S_hat(q,1)-S_hat(q,2),
                    S_hat(q,1),
                    S_hat(q,2)];
            Fel(i,j) = Fel(i,j) - poids(q)*lambda(j)*(invTransBl(2,:)*gradwq(i,:)')*detBl;
        end
    end
end

% A COMPLETER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2023
