% =====================================================
%
% une routine pour la mise en oeuvre des EF P2 Lagrange
% pour l'equation de Laplace suivante, avec conditions de Neumann
%
% | -\Delta u + u= f,   dans \Omega
% |         du/dn = 0,   sur le bord
%
% =====================================================

% lecture du maillage et affichage
% ---------------------------------
nom_maillage = "geomCarreh04.msh";
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri]= lecture_msh_ordre2(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de masse
LL = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for l=1:Nbtri

  % Coordonnees des sommets du triangles
  % A COMPLETER
  S1=Coorneu(Numtri(l,1),:);
  S2=Coorneu(Numtri(l,2),:);
  S3=Coorneu(Numtri(l,3),:);

  % Calcul des matrices elementaires du triangle l
  Kel=matK_elem_p2(S1, S2, S3);
  Mel=matM_elem_p2(S1, S2, S3);

  % On fait l'assemblage des matrices globales
  % A COMPLETER
    for i=1:6 
        for j=1:6
            I=Numtri(l,i);
            J=Numtri(l,j);
            MM(I,J) = MM(I,J) + Mel(i,j);
            KK(I,J) = KK(I,J) + Kel(i,j);
        end
    end 
end

% Calcul du second membre
% -------------------------
% A COMPLETER
% utiliser la routine f.m
FF = zeros(Nbpt,1);
    for i= 1:Nbpt
        FF(i) = f(Coorneu(i,1), Coorneu(i,2));
    end;   
LL = MM*FF;

% Resolution du systeme lineaire
% ----------
UU = (MM+KK)\LL;

% visualisation
% -------------
%affiche_ordre2(UU, Numtri, Coorneu,sprintf('Neumann - %s', nom_maillage));
affiche_ordre2(UU, Numtri, Coorneu,sprintf('Neumann - %s', "f=(1+5*pi^2)*cos(pi*x)*cos(2*pi*y);"));
validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
UUBON = cos(pi*Coorneu(:,1)).*cos(2*pi*Coorneu(:,2))
% Calcul de l'erreur L2
L2 = transpose(UUBON-UU)*MM*(UUBON-UU);
norme2U = transpose(UUBON)*MM*(UUBON);
% Calcul de l'erreur H1
H1 = abs(transpose(UUBON-UU)*KK*(UUBON-UU));
norme1U = abs(transpose(UUBON)*KK*(UUBON));

logerreurL2=log(L2/norme2U);
logerreurH1=log(H1/norme1U);

[logerreurL2,logerreurH1]

h = [log(1/0.4), log(1/0.2), log(1/0.1), log(1/0.05), log(1/0.025)];
yL2 = [-5.1120, -11.7345, -16.6117, -22.1564, -27.2711];
yH1 = [-5.8281, -9.7326, -13.0694, -17.1118, -20.6951];

figure;
loglog(h, yL2, 'b-', 'DisplayName', 'erreur L2'); % Plot de l'erreur L2
hold on;
loglog(h, yH1, 'r-', 'DisplayName', 'erreur H1'); % Plot de l'erreur H1
xlabel("log(1/h)");
ylabel("log(erreurs)");
grid on;
legend('show');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2023
