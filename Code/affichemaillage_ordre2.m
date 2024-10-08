function affichemaillage_ordre2(nom_maillage, titre);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% affichemaillage:
% pour visualiser un maillage triangulaire 2D
%
% SYNOPSIS affichemaillage(nom_maillage, titre)
%
% INPUT  * nom_maillage : le nom d'un fichier de maillage au format msh
%                   AVEC SON SUFFIXE .msh et ENTRE GUILLEMETS (string)
%        * titre (optionel) un titre (string)
%
% OUTPUT une fenetre graphique
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% control on the input args
if (nargin<2), titre = ''; end;

% lecture du fichier nom_maillage.amdba
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri]=lecture_msh_ordre2(nom_maillage); %MODorder2%

%visualisation du maillage
figure;
hold on

% maillage
trimesh(Numtri(:,1:3),Coorneu(:,1),Coorneu(:,2),zeros(Nbpt,1)); %MODorder2%
plot(Coorneu(:,1),Coorneu(:,2),'.r', "markersize", 10); %MODorder2%
view(2);
axis('equal');

% ajouter eventuellement un titre
title(titre);

hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2023
