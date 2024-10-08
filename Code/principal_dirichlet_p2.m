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
nom_maillage = "geomCarreh02.msh";
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
        I=Numtri(l,i);
        for j=1:6
            J=Numtri(l,j);
            MM(I,J) = MM(I,J) + Mel(i,j);
            KK(I,J) = KK(I,J) + Kel(i,j);
        end
    end 
end % for l

% Calcul du second membre
% -------------------------
% A COMPLETER
% utiliser la routine f.m
FF = zeros(Nbpt,1);
    for i= 1:Nbpt
        FF(i) = f(Coorneu(i,1), Coorneu(i,2));
    end;   
LL = MM*FF;
AA = MM+KK;
[tilde_AA, tilde_LL] = elimine(AA, LL, Refneu,Coorneu);

% Resolution du systeme lineaire
% ----------
UU = tilde_AA\tilde_LL;

% visualisation
% -------------
%affiche_ordre2(UU, Numtri, Coorneu,sprintf('Neumann - %s', nom_maillage));
affiche_ordre2(UU, Numtri, Coorneu,sprintf('Dirichlet - %s', "f=(1+2*pi)*sin(pi*x)*sin(pi*y)"));
validation = 'non';
% validation
% ----------
if strcmp(validation,'oui')
UU_exact = sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2))
% Calcul de l'erreur L2
L2 = sqrt(transpose(UU_exact-UU)*MM*(UU_exact-UU))
norme2U = sqrt(transpose(UU_exact)*MM*(UU_exact))
% Calcul de l'erreur H1
H1 = sqrt(abs(transpose(UU_exact-UU)*KK*(UU_exact-UU)))
norme1U = sqrt(abs(transpose(UU_exact)*KK*(UU_exact)))

h=[log10(1/0.4),log10(1/0.2),log10(1/0.1),log10(1/0.05)];
x = h;
yL2=zeros(4);
yH1=zeros(4);
lL2 = [0.0072,5.4134e-04,4.9361e-05,2.7390e-06];
lH1 = [0.0731,0.0112,0.0028,2.6543e-04];
lnorme1U = [4.4350,4.4424,4.4428,4.4429];
lnorme2U = [0.9948,0.9997,1.0000,1.0];
for i=1:4
    yL2(i) = log10(lL2(i)/lnorme2U(i));
    yH1(i) = log10(lH1(i)/lnorme1U(i));
end
loglog(h, yL2, 'b', 'DisplayName', 'L2');
hold on;
loglog(h, yL2, 'b', 'DisplayName', 'L2');
hold on;
loglog(h, yL2, 'b', 'DisplayName', 'L2');
hold on;
loglog(h, yH1, 'r', 'DisplayName', 'H1');
xlabel("log(1/h)");
ylabel("log(rapport des normes)");
grid on;
legend('erreur L2', 'erreur H1');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2023