% TODO
% contrôler que les titres abrégés de chapitres, partie, etc ne dépassent pas 1 ligne
% (faire une commande \widthcheck{maxwidth}{text} )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\PackageName{ouvrage-hermes}
\def\PackageVersion{1.3.pl1}
\def\filedate{2005/11/14}

%%% Prélude %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=== Identification Part

\NeedsTeXFormat{LaTeX2e}

\def\fileversion{\PackageVersion}
\ProvidesClass{\PackageName}[
     \filedate \space Version:\space \PackageVersion \space ]
\typeout{ouvrage-hermes.cls: Latex style for Hermes books and treatises, Roger Rousseau, 2005}
\ProvidesFile{\PackageName.cls}[\filedate \space Version:\space \PackageVersion]

%=== Initial Code Part

% Définition des variables d'option et leur valeur par défaut
% Les tests if ne comportent pas de `@' pour pouvoir être utilisé dans les documents
% En particulier \iftreatise peut être utile

\newif\iftreatise     \treatisefalse     % Cas particulier des traités
\newif\iffrenchlang   \frenchlangtrue    % ouvrage en français, false => in english
\newif\ifcropmarks    \cropmarksfalse    % affichage du contour de page et de corps du texte
\newif\ifallpages     \allpagesfalse     % imprime aussi les pages blanches
\newif\ifemptyheaders \emptyheadersfalse % aucun entête de page, mais numérotation hors cadre

%=== Declaration of options

\DeclareOption{treatise}{\treatisetrue}
\DeclareOption{english}{\frenchlangfalse}
\DeclareOption{cropmarks}{\cropmarkstrue}
\DeclareOption{allpages}{\allpagestrue}
\DeclareOption{empty}{\emptyheaderstrue}

\DeclareOption*{\typeout{*** Unknown option ``\CurrentOption'' ****\\
Syntax: documentclass[treatise,english,cropmarks,allpages,empty]{ouvrage-hermes}}}
\PassOptionsToClass{10pt,twoside,fleqn}{book}

%=== Execution of options

\ExecuteOptions{10pt}         % Defaults options
\ProcessOptions               % Évaluation des déclarations d'options

\typeout{*** current options for ouvrage-hermes.cls ************}

\iftreatise     \typeout{type of work  : treatise option}%
\else           \typeout{type of work  : book (by default)}%
\fi

\iffrenchlang   \typeout{main language : french (by default)}%
\else           \typeout{main language : english option}%
\fi

\ifcropmarks    \typeout{page style    : cropmarks option}%
\else           \typeout{page style    : no cropmarks (by default)}%
\fi

\ifemptyheaders \typeout{headers       : no (empty option)}%
\else           \typeout{headers       : yes (by default)}%
\fi

\ifallpages     \typeout{white pages   :  printing (allpages option)}%
\else           \typeout{white pages   : no printing (by default)}%
\fi

\typeout{*******************************************************}

%=== Package loading

\LoadClass{book}
% Attention ! chapterbib doit être chargé avant babel !
\iftreatise
\RequirePackage{chapterbib}
\fi

% On n'utilise plus french de B. Gaulle, mais frenchb (de babel)
\RequirePackage[francais,english]{babel}

% Attention ! avec babel, tout appel à \selectlanguage active un espacement spécial pour
% les listes : nous ne nous laisserons pas faire !
\let\ORIselectlanguage\selectlanguage
\gdef\selectlanguage#1{\ORIselectlanguage{#1}\@mkhermes}

\RequirePackage{ifthen} % pas indispendable, mais standard
\RequirePackage{xspace}
  % xspace est utilisé par babel pour améliorer l'utilisation de l'espace après les macros
  % comme \no, 1\up{er}...

\RequirePackage{makeidx} 
\RequirePackage{multicol} % pour l'index (\twocolumn ne traite pas bien le \enlargethispage)
\RequirePackage{epsfig}   % pour la commande \EpsFigure
\RequirePackage{textcomp}

%=== Packages Options as variables

% Les options de Babel sont donc appelées dans \AtBeginingDocument.

%%% Fin du prélude

%%% Code des Commandes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Changement des options en cours d'ouvrage (très spécial)

\newcommand{\CropMarksOn}{\gdef\@DrawPageFrame{\@PageFrame}}
\newcommand{\CropMarksOff}{\gdef\@DrawPageFrame{\hspace*{0pt}}}

%%% Primitives internes ----------------------------------------------------------------

\def\@BS{\symbol{"5C}}
\def\@ARGS{\{\ldots\}}
\def\@ARG#1{\{\normalfont\textit{#1}\}}

%%% FONTES : de la forme FonteNom

\def\FonteTexte{%
   \normalfont\fontsize{10}{12}\selectfont} % Police du texte normal

\def\FonteEnTete{%
   \normalfont\fontsize{9}{12}\selectfont}% Police des en-têtes

\def\FonteSousTitreChapitre{%
   \normalfont\fontsize{14}{22}\selectfont}

\def\FonteTitreChapitre{%
   \normalfont\fontsize{18}{22}\selectfont}

\def\FonteNumeroPartie{%
   \rmfamily\mdseries\upshape\fontsize{14pt}{22pt}\selectfont}%

\def\FonteTitrePartie{%
   \rmfamily\mdseries\upshape\fontsize{18pt}{22pt}\selectfont}%

\def\FonteSectionI{%
   \normalfont\normalsize\bfseries}%

\def\FonteSectionII{%
   \normalsize\normalfont\bfseries\itshape}%

\def\FonteSectionIII{%
   \normalsize\normalfont\itshape}%

\def\FonteSectionIV{%
   \normalsize\normalfont}%

% On laisse footnotesize à sa valeur usuelle de normalsize - 2 pts, soit 8 pts.
% Ainsi FonteNoteBasPage sera la référence pour 9/11pts
% L'appel de note doit être dans une taille relative à la police courante et non
% dans l'absolu à 8 pts
\def\FonteAppelNote{\protect\normalfont\footnotesize}

\def\FonteNoteBasPage{%
   \normalfont\fontsize{9pt}{11pt}\selectfont}

\def\FonteTypeLegende{% Par exemple Tableau, Figure, Programme...
   \FonteNoteBasPage\bfseries}

\def\FonteTitreLegende{%
   \FonteNoteBasPage\itshape}

\def\FonteBibliographie{\rmfamily\mdseries\upshape\fontsize{9}{11}\selectfont}

%%% VARIABLES : de la forme VAR@Nom

\def\@VARTitreComplet{} 
% Titre complet de l'ouvrage tel qu'il se présente sur la couverture
% identique à \@title, mais ce dernier est redéfini dans les macros de latex ou book

\def\VAR@TitreAbrege{%
% Titre abrégé de l'ouvrage pour les en-tête de pages de gauche (paires)
  Titre de l'ouvrage, à définir par \texttt{\@BS title[titre abrégé]\{titre\} }}

\def\@VARAuteurs{}
% Auteurs du livre ou coordinateur d'un ouvrage collectif ou éditeur scientifique...
% Identique à \@author, mais cette dernière est effacée par latex ou book après usage

\def\@author{%
% Auteurs du livre ou coordinateur d'un ouvrage collectif ou éditeur scientifique...
   Auteur(s)\\
   {\FonteTexte\textit{A définir par la commande \texttt{\upshape\@BS author\@ARGS}.
    Dans le cas de plusieurs auteurs, séparer chaque auteur par 
    \texttt{\upshape\@BS andauthor}.
}}}

\def\VAR@ChapitreAbrege{%
  \FonteEnTete%
  Titre abrégé du chapitre, à définir par \texttt{\@BS chapter[titre abrégé]\{titre\}}}

\def\VAR@Phone{} % Des Auteurs pour le service de fabrication...
\def\VAR@Fax{}
\def\VAR@Email{}

% Présentation --------------------------------------------------------------------------

% Polices du corps du texte

% L'encodage recommandé est T1 et saisie du texte directement en latin1
  \usepackage[T1]{fontenc}
  \usepackage[latin1]{inputenc}%C. QUEINNEC

% Latex/Unix: Times romain 
\renewcommand{\rmdefault}{ptm}  % Postscript Times New Roman
\renewcommand{\sfdefault}{phv}  % Postscript Helvetica
%\renewcommand{\ttdefault}{pcr} % Police Courier : un peu large
\renewcommand{\ttdefault}{cmtt} % Police Latex tt : plus jolie et étroite que Courier
                                 % mais n'existe pas dans toutes les tailles
\renewcommand{\bfdefault}{b}    % le bold extra (bx) n'existe pas en Times

\newlength{\interligne}       % Ne change pas dans ce style: quand on parle
\setlength{\interligne}{12pt} % d'une ligne blanche, c'est ça : 4.22mm
% On peut utiliser cette commande dans le texte, par exemple
% \setlength{\baselineskip}{1\interligne} pour revenir à l'espacement normal
\def\@espace{\hspace*{0.5em}}  % Hermès parle parfois d'espaces (une espace)
                                        % pour les écarts horizontaux ~= 1/2 em

% Format des pages ----------------------------------------------------------------------

% Format des pages (pour papier a4 ou letter)
% Rappel : A4 = 210mm x 297 mm USLETTER = 8.5in x 11in = 215.9mm * 279.4mm
% soit USLETTER = 6mm de plus dans la marge droite et 17.6mm de moins dans le bas de page

% 1pt = 1/72,27 inch= 25,4/72.27mm = 0.35146 mm
% 1mm = 2.84528 pt

\setlength{\textwidth}{120mm}  % 120mm = 341pt
\setlength{\textheight}{183mm} % 183mm = 521pt (43.4 lignes * 12pt)

%--- Dimensions horizontales

\setlength{\hoffset}{-1in}            % à ajuster par dvips selon imprimante
\setlength{\oddsidemargin}{45mm}      % pour que la marge de gauche
\setlength{\evensidemargin}{45mm}     % soit exactement de 45mm
  % Largeur de la marge de droite pour papier:
  %  - A4     : largeur papier = 210 mm => 210 - (120+45)  = 45  mm
  %  - Letter : 8.5in = 8.5*25.4mm = 215.9mm => 215.9-165  = 50.9mm

%--- Dimensions verticales

\setlength{\voffset}{-1in}       % Élimination de la marge par défaut de 1 pouce

\setlength{\topmargin}{47mm}     % distance bord sup feuille au haut de l'en-tête
\settoheight{\headheight}{2.5mm} % hauteur de l'en-tête
\setlength{\headsep}{7.5mm}      % marge entre corps du texte et en-tête
                                 % => 57mm = 47mm + 2.5mm + 7.5mm = marge au dessus du
                                 % texte d'une feuille A4 non massicotée

%---------------------------------------------------------------------------------------------
% Tracé du cadre indiquant les bords du papier une fois l'ouvrage massicoté
% au format des ouvrages Hermès environ (160mm x 240mm) ou (155mm x 235mm)

% Ce tracé est obtenu par l'option \documentclass[cropmarks]{ouvrage-hermes}
% (par défaut pas de tracé) et est commandé à partir de l'en-tête (picture).
% L'avantage de ce système est de ne pas influencer la mise en page de Latex

%-- Primitives internes

\def\@DrawPageFrame{\hspace*{0pt}}

\def\@PageFrame{%
   % Tracé des contours de texte et de page pour la mise au point des problèmes
   % de débordement indiqués par LaTeX
   \thinlines%
   \setlength{\unitlength}{1mm}%
   \begin{picture}(0,0)%
   % Les mesures sont faites à partir du point de référence qui est posé à gauche
   % de l'en-tête (7.5mm au dessus du corps du texte)
   % L'ouvrage massicoté est de 160mmx240mm (ou 155mm x 235mm)
   % Le corps du texte est inscrit dans un rectangle de 120mmx183mm
   % Les marges du corps du texte dans une feuille A4 sont de :
   %  - haute et basse de 57 mm (297-183)/2
   %  - gauche et droite de 45 mm (210-120)/2

%%%% Cadre des bords du papier massicoté : cas de pages en 155x235
% En format A4, les bandes supprimées au dessus sont de (297-235)/2=31mm
% donc la distance du coin supérieur gauche au dessus du point de référence est égale à
%    57-(31+7.5) = 18,5mm = distance entre point de référence et bord supérieur papier
   \put(-17.5,18.5){\line(1,0){155}}    % 18.5mm = 57 - (31 + 7.5) %
                                        % 17.5mm = (155 - 120)/2   %
   \put(-17.5,-216.5){\line(1,0){155}}  % 216.5  = 235-18.5        %
   \put(-17.5,18.5){\line(0,-1){235}}   % 235                      %
   \put(137.5,18.5){\line(0,-1){235}}   % 137.5  = 120 + 17.5mm    %

%%%% Cadre des bords du papier massicoté : cas de pages en 160x240
% En format A4, les bandes supprimées au dessus sont de (297-240)/2=28.5mm
% donc la distance du coin supérieur gauche au dessus du point de référence est égale à
%    57-(28.5+7.5) = 21mm = distance entre point de référence et bord supérieur papier
%   \put(-20,21){\line(1,0){160}}      % 21mm = 57 - (28.5 + 7.5) %
%                                      % 20  mm = (160 - 120)/2   %
%   \put(-20,-219){\line(1,0){160}}    % 219  = 240-21            %
%   \put(-20,21){\line(0,-1){240}}     % 240                      %
%   \put(140,21){\line(0,-1){240}}   % 140  = 120 + 20mm      %
%%%% Cadre des bords du texte
   \put(0,-7.5){\line(1,0){120}}      % marque du bord supérieur du texte
   \put(0,-190.5){\line(1,0){120}}    % 7.5+183 = 190.5mm: bas du texte
   \put(0,-7.7){\line(0,-1){183}}     % côté gauche du texte
   \put(120,-7.5){\line(0,-1){183}}   % côté droit du texte
   \end{picture}}%

\ifcropmarks \CropMarksOn \else \CropMarksOff \fi
\flushbottom

%--- Définition des styles de page (du point de vue des en-têtes et pieds de page)
%
% Critique : On pourrait utiliser les noms classiques plain et empty
%%%
% Pas de pied de page pour les ouvrages: 
% mais on l'utilise pour numéroter les pages en sortant du cadre massicoté

% On redéfinit \clearpage et \doubleclearpage pour faire apparaître
% des pages blanches sans en-tête (sauf numéro de page hors limites)
% mais avec le cadre si l'option cropmarks est demandée

\def\clearpage{% latex.ltx + \thispagestyle
  \ifvmode                          %
    \ifnum \@dbltopnum =\m@ne       %
      \ifdim \pagetotal <\topskip   %
        \hbox{}                     %
      \fi\fi\fi                     %
  \newpage\thispagestyle{emptypage} %
  \SmallFnRule                      %
  \write\m@ne{}                     %
  \vbox{}                           %
  \penalty -\@Mi}

\def\cleardoublepage{% latex.ltx + 2 modifs :
    % \thispagestyle
    % 2e page blanche seulement si l'option allpages est armée.
    \clearpage\if@twoside \ifodd\c@page\else %
    \ifallpages %
    \hbox{}\newpage\thispagestyle{emptypage}\if@twocolumn\hbox{}\newpage\fi%
    \else \addtocounter{page}{1}\fi\fi\fi}

\let\ORImaketitle\maketitle
\def\maketitle{\ORImaketitle\cleardoublepage}
 
% Presque \ps@empty

\def\ps@emptypage{%
   % Une première page est normalement impaire, mais on définit les deux...
   \renewcommand{\@oddhead}{\parbox{\textwidth}{\@DrawPageFrame}}%
   \renewcommand{\@evenhead}{\parbox{\textwidth}{\@DrawPageFrame}}%
   \setlength{\footskip}{40mm} % réduire à 33mm pour USLETTER
   \renewcommand{\@evenfoot}{\FonteEnTete \hfill \thepage\hfill}%
   \renewcommand{\@oddfoot}{\FonteEnTete \hfill \thepage\hfill}}% emptypage style

% \ps@empty

\def\ps@empty{%
   \renewcommand{\@oddhead}{}%
   \renewcommand{\@evenhead}{}%
   \setlength{\footskip}{40mm}%
   \renewcommand{\@evenfoot}{}%
   \renewcommand{\@oddfoot}{}}

% Équivalent de \ps@plain

\newlength{\@ESPage}\setlength{\@ESPage}{4mm}% 5 espaces blancs pour entete

\def\ps@plainpage{% Page normale avec en-têtes (alias normalpage)
   \ifemptyheaders % pas d'en-tête
      \renewcommand{\@evenhead}{\parbox{\textwidth}{\@DrawPageFrame}}%
      \renewcommand{\@oddhead}{\parbox{\textwidth}{\@DrawPageFrame}}%
      \setlength{\footskip}{40mm} % réduire à 33mm pour USLETTER
      \renewcommand{\@evenfoot}{\FonteEnTete \hfill \thepage \hfill}%
      \renewcommand{\@oddfoot}{\FonteEnTete \hfill \thepage \hfill}%
   \else % En-tête normal
      \ifemptyheaders\else % Impression En-tête
         \renewcommand{\@evenhead}{%
            \parbox{\textwidth}{%
                \@DrawPageFrame\FonteEnTete\thepage\hspace{\@ESPage}\VAR@TitreAbrege}}%
         \renewcommand{\@oddhead}{%
            \parbox{\textwidth}{%
                \@DrawPageFrame\hspace*{\fill}\FonteEnTete\VAR@ChapitreAbrege%
                \hspace{\@ESPage}\thepage}}%
         \renewcommand{\@evenfoot}{}%
         \renewcommand{\@oddfoot}{}%
        \fi %
   \fi} % plainpage style
\pagestyle{plainpage}% doit être placé après les définitions précédentes

% Titre de l'ouvrage --------------------------------------------------------------------

\renewcommand{\title}[2][]{
% usage: \title[shorttitle]{title}
   \ifthenelse{\equal{#1}{}}{%Y
      \@shorttitle{#2}%
   }{%
      \@shorttitle{#1}%
   }
   \gdef\@title{#2}
   \gdef\@VARTitreComplet{#2}
}

%--- titre courant (ou abrégé) pour les en-têtes de pages de droite (paires)

\newsavebox{\@shtitle}
\def\@shorttitle#1{%
   % #1: Titre abrégé de l'ouvrage (moins de 40 caractères, pour en-tête)
   % test de la largeur < 55mm
   \sbox{\@shtitle}{\FonteEnTete #1}%
    \ifthenelse{\lengthtest{\wd\@shtitle > 55mm}}{%
     \@latex@warning{short title defined by \@backslashchar title[short]{long} too long:
       (55mm max)}
    }{}
    \gdef\VAR@TitreAbrege{#1}
}

%--- Titre de l'ouvrage

\def\@Titre#1{{%
   % #1: Titre complet de l'ouvrage (peut incorporer des `\\')
   \parskip \z@%
   \FonteTitre #1%
   \par\vspace*{20.4pt}% 2 lignes = 24 pts - (21.6-18)
}}

%--- Auteurs : distinguer Auteurs d'un livre et Auteurs d'un chapitre

\def\author#1{\gdef\@author{#1}\gdef\@VARAuteurs{#1}}

%------ Commande de sectionnement : \frontmatter,\mainmatter,\backmatter -----------------

\renewcommand\frontmatter{% adaptation de book.cls
  % Sert aussi pour initialiser des variables aprés \begin{document}
  % Utilisation obligatoire dans tout ouvrage
  \cleardoublepage %
  \@mainmatterfalse %
  \pagenumbering{arabic} %
}

\renewcommand\mainmatter{% adaptation de book.cls
    \cleardoublepage
  \@mainmattertrue
  %  \pagenumbering{arabic} % Remet à 0 le numéro de page
}

\renewcommand\backmatter{%% adaptation de book.cls
    \cleardoublepage %
    \@mainmatterfalse %
    \@endofappendices %
}

%------ Espaces verticaux pour les pages spéciales
% 1pt = 25.4/72.27 mm = 0.35146 mm (point typographique)
% 60 pt = 21.1 mm
% 90 pt = 31.6 mm
% 17 pt =  6   mm
% On utilise des espaces différents pour corriger les décalages parasites
% mais le résultat est parfait et stable

\newlength{\@DessusTitrePT}\setlength{\@DessusTitrePT}{11.2mm}
\newlength{\@DessusTitrePTS}\setlength{\@DessusTitrePTS}{21.8mm} %9.2mm+22pt+17pt
\newlength{\@DessusTitreCN}\setlength{\@DessusTitreCN}{49pt} % 60pt - 11pt
\newlength{\@DessusTitreCS}\setlength{\@DessusTitreCS}{40pt} % 60pt - 20pt

\newlength{\@DessousTitreCN}\setlength{\@DessousTitreCN}{79pt} % 90pt - 11pt
\newlength{\@DessousTitreCS}\setlength{\@DessousTitreCS}{79pt} % 90pt - 19pt

\newlength{\@Intertitre}         % Esp. vert entre un type de page et un titre
\setlength{\@Intertitre}{17pt}   % = 6mm

%%%% ======== Part (Niveau -1) ====================================================

\renewcommand{\part}{% adaptation de book.cls
  \setcounter{footnote}{0} % utile seulement pour les pages exceptionnelles avec footnote
  \if@openright\cleardoublepage\else\clearpage\fi%
  \thispagestyle{emptypage}    %
  \enlargethispage{-3pt}       % correction d'un décalage inexpliqué du texte sur page titre
  \if@twocolumn                %
    \onecolumn                 %
    \@tempswatrue              %
  \else                        %
    \@tempswafalse             %
  \fi                          %
  \secdef\@part\@spart         %
}

\def\@part[#1]#2{% adaptation de book.cls
   % #1 : Titre abrégé de la partie pour la table des matières
   % #2 : Titre complet de la partie
   %
   % Mise à jour du compteur de partie et de de la table des matières
    \vspace*{\@DessusTitrePT}           %
    \ifnum \c@secnumdepth >-2\relax     %
      \refstepcounter{part}             %
      \addcontentsline{toc}{part}{\textsc{\OrdinalInWords{part}\space\@partname.~#2}}%
    \else                               %
      \addcontentsline{toc}{part}{\textsc{#2}} %
    \fi                                 %
      \vspace*{0pt}                     %
      \begin{flushright}                %
       \interlinepenalty \@M            %
       \ifnum \c@secnumdepth >-2\relax  %
         {\FonteNumeroPartie\textsc{\OrdinalInWords{part}~\@partname}}\\[\@Intertitre] %
       \fi                              %
       {\FonteTitrePartie #2}           %
     \end{flushright}                   %
     \@endpart}                         %

\def\@spart#1{% adaptation de book.cls
   % Affichage d'une partie sans numérotation (\part*)
    \vspace*{\@DessusTitrePTS}          %
    \begin{flushright}                  %
     \interlinepenalty \@M              %
     {\FonteTitrePartie #1}             %
    \end{flushright}                    %
     \@endpart}                         %

\def\@endpart{% adaptation de book.cls
   % Terminaison de la page de partie après affichage
   \cleardoublepage                     %
   \if@tempswa                          %
     \twocolumn                         %
   \fi}                                 %

%%%% ========== Chapter (Niveau 0) =======================================================

\renewcommand{\chapter}{% adaptation de book.cls
  \if@openright\cleardoublepage\else\clearpage\fi %
  \thispagestyle{emptypage}%                      %
  \enlargethispage{-3pt}       % correction d'un décalage inexpliqué du texte sur page titre
  \global\@topnum\z@                              %
  \@afterindentfalse                              %
  \secdef\@chapter\@schapter}                     %

\def\@chapter[#1]#2{% adaptation de book.cls
% #1 : Titre abrégé de la partie pour la table des matières
% #2 : Titre complet de la partie
%
% Mise à jour du compteur de chapitre et de de la table des matières
    \ifnum \c@secnumdepth >\m@ne                  %
       \if@mainmatter % chapitre numéroté de la partie centrale
          \refstepcounter{chapter}                %
          \typeout{\@chapapp\space\thechapter.}   %
          \addcontentsline{toc}{chapter}{\@chapapp\space\thechapter.~#2}%
       \else % chapitre non numéroté, mais dans TOC, comme Introduction, Bibliographie
             % situé soit dans \frontmater, soit dans backmatter
          \addcontentsline{toc}{chapter}{#2}%
          \vspace*{-9pt} % correction nécéssaire
       \fi %
    \else  % demande explicite de non numérotation par la valeur de \secnumdepth
          \addcontentsline{toc}{chapter}{#2}%
    \fi %
    \gdef\VAR@ChapitreAbrege{#1}  %
    \addtocontents{lof}{\protect\addvspace{11\p@}}%
    \addtocontents{lot}{\protect\addvspace{11\p@}}%
    %@@ \addtocontents{lop}{\protect\addvspace{11\p@}}% liste des programmes
    % Affichage du titre et noms d'auteurs éventuels
    \ifthenelse{\equal{\AuthorOfChapter}{}}{}{%else %
          \def\and{, }% on redéfinit \and à `,'
          \addtocontents{toc}{\protect\TextInToc{0}{0pt}{\AuthorOfChapter}}%
          \@defand% on restaure la définition de \and selon la langue (et, and...)
        }%
    \addtocontents{toc}{\protect\addvspace{5pt plus 1pt minus 1pt}}% avant 1ère section
    \if@twocolumn %
        \@topnewpage[\@makechapterhead{#2}]%
    \else %
        \@makechapterhead{#2}%
        \@afterheading% traitement du premier paragraphe
    \fi}%

\def\@makechapterhead#1{% adaptation de book.cls
% Affichage du titre du chapitre numéroté et ne figurant ni dans TOC, ni dans les En-têtes
% de page. Exemples : préface, table des matières (mais avec en-tête)...
  \vspace*{\@DessusTitreCN}           %
  {\parindent \z@ \iftreatise\centering\else\raggedleft\fi \normalfont \parskip \z@ %
    \ifnum \c@secnumdepth >\m@ne      %
      \if@mainmatter                  %
        {\FonteSousTitreChapitre \@chapapp\space\thechapter}\\[\@Intertitre] %
      \fi                                   %
    \fi                                     %
    \interlinepenalty\@M                    %
    {                                       %
    \if@mainmatter                          %
        \FonteTitreChapitre                 %
    \else                                   %
        \FonteSousTitreChapitre             %
    \fi                                     %
    #1\\[\@DessousTitreCN]}%
}}%

\def\@makeschapterhead#1{% adaptation de book.cls           %
% Affichage du titre de \chapter*                           %
  \gdef\VAR@ChapitreAbrege{#1}                    %
  \vspace*{\@DessusTitreCS}                                 %
  {\parindent \z@ \iftreatise\centering\else\raggedleft\fi \normalfont \parskip \z@\relax %
    \interlinepenalty\@M                                    %
    {\FonteSousTitreChapitre                                %
    \VAR@ChapitreAbrege\\[\@DessousTitreCS]}}}              %

%%%% ========== ChapterAuthor (Niveau 0) =================================================

% L'utilisation de \bsc est recommandée, mais doit être utilisée avec précaution pour les
% noms composés : risque de pbs de justification
\newcommand{\Name}[1]{\textsc{#1}}

\def\AuthorOfChapter{}
\def\@PresentationOfAuthors{}
\newcommand{\PresentationOfAuthors}[1]{\gdef\@PresentationOfAuthors{#1}}

\newcommand{\ChapterAuthor}[3][]{%
% #1 : Titre abrégé du chapitre pour la table des matières
% #2 : Titre complet du chapitre
% #3 : Noms des Auteurs du chapitre, selon le format :
%      Prénom1 NOM1, Prénom2 NOM2, ...\and Prénom1 NOM1
  \gdef\AuthorOfChapter{#3}
  \ifthenelse{\equal{#1}{}}{\chapter[#2]{#2}}{\chapter[#1]{#2}}
  \begingroup
     \ifthenelse{\equal{\@PresentationOfAuthors}{}}%
      {%Présentation par défaut
        \iffrenchlang
        \@sfootnotetext{Chapitre rédigé par #3.}%
        \else
        \@sfootnotetext{Chapter written by #3.}%
        \fi
      }{% Sinon Présentation fournie
        \@sfootnotetext{\@PresentationOfAuthors #3.}
      }
  \endgroup%
  \gdef\AuthorOfChapter{}%
  \PresentationOfAuthors{}%
}

% Abréviations
\let\Chapter\ChapterAuthor
\let\Author\AuthorOfChapter

%%%% ========== Annexes (Niveau 0) =================================================

\newif\if@appendix \@appendixfalse
\let\Orichapter\chapter

\def\@endofappendices{
   % Pour signaler la fin des annexes.
   % utilisée par \backmatter, \thebibliography, \index
   \@appendixfalse
   \let\chapter\Orichapter
}


\renewcommand\appendix{% adaptation de book.cls
  \@appendixtrue
  \par
  \setcounter{chapter}{0}%
  \setcounter{section}{0}%
  % Inutile de modifier tocdepth : cela ne sert que dans \tableofcontents placé
  \gdef\@chapapp{\appendixname}%
  \gdef\thechapter{\@Alph\c@chapter}
  \addcontentsline{toc}{chapter}{\appendicesname}%
  \addtocontents{toc}{\protect\vspace{5pt}}% avant les annexes
  \gdef\chapter{%
     \if@openright\cleardoublepage\else\clearpage\fi
     \thispagestyle{emptypage}%
     \enlargethispage{-3pt}       % correction d'un décalage inexpliqué du texte sur page titre
     \global\@topnum\z@
     \@afterindentfalse
     \secdef\@Appendix\@sAppendix}%
}
\def\@Appendix[#1]#2{%
   \refstepcounter{chapter}%
   \addcontentsline{toc}{appendix}%
      {\protect\numberline{\thechapter.}~#2}%
   \gdef\VAR@ChapitreAbrege{#1}
    \if@twocolumn
        \@topnewpage[\@makechapterhead{#2}]%
    \else
        \@makechapterhead{#2}%
        \@afterheading% traitement du premier paragraphe
    \fi
}

\def\@sAppendix#1{%
 \if@twocolumn
   \@topnewpage[\@makeschapterhead{#1}]%
   \else
   \@makeschapterhead{#1}%
   \@afterheading
   \fi}

%--- Traitement des sections dans les annexes

% Les sections des annexes sont présentées comme dans la partie centrale, mais
% apparaîssent différemment dans la table des matières.
% Je redéfinis \@sect de latex.ltx pour détecter l'apparition précédente de \appendix
% et appeler des styles de lignes de toc différentes selon le cas.
% C'est de loin la solution la plus simple et la plus efficace.

\def\@sect#1#2#3#4#5#6[#7]#8{% adaptation de la version de latex.ltx
% Il s'agit de traiter différemment le placement dans toc selon que l'on a franchi ou non
% la commande \appendix. 
% J'ajoute un point et un espace après le numéro de section dans la toc.
% Le reste est inchangé, car il est difficile de traiter
% correctement en latex les espacements verticaux lorsque deux sections se suivent...
  \ifnum #2>\c@secnumdepth
    \let\@svsec\emptyheaders
  \else
    \refstepcounter{#1}%
    \protected@edef\@svsec{\@seccntformat{#1}\relax}%
  \fi
  \@tempskipa #5\relax
  \ifdim \@tempskipa>\z@
    \begingroup
      #6{%
        \@hangfrom{\hskip #3\relax\@svsec}%
          \interlinepenalty \@M #8\@@par}%
    \endgroup
    \csname #1mark\endcsname{#7}%
% Deux cas différents selon \appendix ou non
   \if@appendix
      \ifnum #2=1 % format particulier pour les sections en annexe
         \addcontentsline{toc}{appendixi}{%
             \ifnum #2>\c@secnumdepth \else
             \protect\numberline{\csname the#1\endcsname.}~%
             \fi
          #7}%
      \fi
      \ifnum #2=2 % format particulier pour les sous-sections en annexe
         \addcontentsline{toc}{appendixii}{%
             \ifnum #2>\c@secnumdepth \else
             \protect\numberline{\csname the#1\endcsname.}~%
             \fi
          #7}%
      \fi
      % Pour les autres sections en annexe, les chapitres sont traités dans \
      % et les autres sections (niveaux 3, 4, 5) ne doivent pas être cités dans toc
   \else % traitement normal de latex.ltx
       \addcontentsline{toc}{#1}{%
          \ifnum #2>\c@secnumdepth \else
%            \protect\numberline{\csname the#1\endcsname}%
            \protect\numberline{\csname the#1\endcsname.}~%
          \fi
          #7}%
   \fi
  \else
    \def\@svsechd{%
      #6{\hskip #3\relax
      \@svsec #8}%
      \csname #1mark\endcsname{#7}%
      \addcontentsline{toc}{#1}{%
        \ifnum #2>\c@secnumdepth \else
%          \protect\numberline{\csname the#1\endcsname}%
          \protect\numberline{\csname the#1\endcsname.}~%
        \fi
        #7}}%
  \fi
  \@xsect{#5}}

% Intertitres (sections) niveau 1 à 4 -----------------------------------------------

% cf \Latex Companion 2.3 pp 18-31
% Les définitions suivantes visent à donner la police particulière
% exigée par Hermes uniquement dans les titres de section et non pas dans les
% référence. Ainsi, 2.2. apparaîtra en gras droit dans le titre et en
% police normale dans les références par \label ... \ref.
% Si l'on veut la même présentation dans les références de sections,
% il faut tout simplement incorporer la définition voulue dans celle
% de \thesection ... \thesubsection

%%%----------------------------------------------------------------------------------------

% J'utilise \@startsection pour les définitions des sections des niveaux 1 à 4. (Latex Companion)
% Attention, il faut tenir compte de \parskip dans les dimensions de \beforeskip et \afterskip.

% J'ai dû modifier \@sect de latex.ltx pour ajouter un point après le numéro de
% dans la table des matières, car les numéros comme \thesection, \thesubsection ne doivent
% pas comporter de point final pour les \ref{...} : voir redéfinition plus haut pour les annexes.

% J'ai dû redéfinir \@xsect de latex.ltx qui traite le afterskip, pour pouvoir
% traiter correctement le afterskip négatif de -6pt, pour avoir un afterskip effectif de 5 pt
% puisque \parskip est de 11pt.


\def\@xsect#1{% Adaptation de latex.ltx : traitement de afterskip
  % Accepte un afterskip négatif comme déplacement vertical
  % Mais n'interprète plus les afterskip négatifs  comme une demande de suppression 
  % du passage à la ligne (run-in heading).
    \par \nobreak
    \vskip #1
    \vspace*{-11pt} % Bricolage ! Attention 11pt=\parskip
    \@afterheading
    \ignorespaces%
}

%%%-----------------------------------------------------------------------------------------

% On numérote jusqu'à la section de niveau 4
\setcounter{secnumdepth}{4}

% Pour les titres de section sur plusieurs lignes, Hermès ne veut pas d'indentation
\def\@hangfrom#1{%
       \setbox\@tempboxa\hbox{{#1}}%
%       \hangindent \wd\@tempboxa\noindent\box\@tempboxa %% ORIGINAL de latex.ltx
       \hangindent \wd 0 \noindent\box\@tempboxa %% MODIF
}

% Attention : les dimensions des interlignes before et after skip
% doivent tenir compte du parskip qui est de 6pts et qui s'ajoute aux valeurs
% indiquées

%--- Niveau 1: section : gras,  beforeskip = 2 lignes, afterskip = 1 ligne

\renewcommand{\thesection}{\thechapter.\arabic{section}}%

\renewcommand{\section}{%
   \renewcommand{\@seccntformat}[1]{\thesection.\hspace{0.5em}}%Doit être défini ici, en interne
   \@startsection{section}%  % nom de l'inter
   {1}%                      % niveau de l'inter
   {0pt}%                    % l'indentation du titre et du texte suivant
   {11pt plus 3pt minus 2pt} % beforeskip : +parskip = 22pt
%   {10pt plus 2 pt minus 1pt}  % afterskip = 11pt.
   {11pt}                    % afterskip = 11pt.
   {\FonteSectionI}
}

%--- Niveau 2: sous-section : gras+ital, beforeskip = 2 lignes, afterskip = 1 ligne

\renewcommand{\thesubsection}{\thesection.\arabic{subsection}}%
\renewcommand{\subsection}{%
  \renewcommand{\@seccntformat}[1]%
               {{\FonteSectionI\thesubsection.\hspace{0.5em}}}%
  \@startsection%
   {subsection}%            % nom de l'inter
   {2}%                     % niveau de l'inter
   {0pt}%                   % l'indentation du titre et du texte suivant
   {11pt plus 2pt minus 2pt}  % beforeskip : +parskip = 22pt
   {11pt}                     % afterskip = 11pt.
   {\FonteSectionII}}%      % style

%--- Niveau 3: sous-sous-section : italique, beforeskip=1 lignes, afterskip=1 ligne

\renewcommand{\thesubsubsection}{\thesubsection.\arabic{subsubsection}}%
\renewcommand{\subsubsection}{%
  \renewcommand{\@seccntformat}[1]{%
               {\normalsize\normalfont\thesubsubsection.\hspace{0.5em}}}%
  \@startsection%
   {subsubsection}%         % nom de l'inter
   {3}%                     % niveau de l'inter
   {0pt}%                   % l'indentation du titre et du texte suivant
   {0pt plus 1pt minus 1pt} % beforeskip  + parskip = 11pt.
   {5pt}%                   % afterskip = 5pt.
   {\FonteSectionIII}}  % style

%--- Niveau 4: paragraphes: police normale, beforeskip=1 ligne, afterskip=1 ligne

\renewcommand{\theparagraph}{\thesubsubsection.\arabic{paragraph}}%
\renewcommand{\paragraph}{%
  \renewcommand{\@seccntformat}[1]{%
               {\normalsize\normalfont\theparagraph.\hspace{0.5em}}}%
  \@startsection%
   {paragraph}%             % nom de l'inter
   {4}%                     % niveau de l'inter
   {0mm}%                   % l'indentation du titre et du texte suivant
   {0pt plus 1pt minus 1pt} % beforeskip  + parskip = 11pt.
   {5pt}%                   % afterskip = 5pt.
   {\FonteSectionIV}        %
} % style

%--- Niveau 5: interdit
\renewcommand{\subparagraph}{%
    \@latex@error{subparagraphs forbidden in hermes-ouvrage.cls}\@ehd}

% Paragraphes (ou alinéas) ----------------------------------------------------------------

\setlength{\parindent}{5mm}                      % renfoncement (ou retrait ou alinea positif)
\newlength{\@normalparskip}\setlength{\@normalparskip}{11pt plus 3pt minus 2pt}
 % espacement entre les paragraphes
\setlength{\parskip}{\@normalparskip}

% Disable single lines at the start of a paragraph, at the end of a page
% (lignes orphelines)
\clubpenalty = 10000
% Disable single lines at the end of a paragraph, at the beginning of a page
% (lignes veuves)
\widowpenalty = 10000
%\displaywidowpenalty = 10000
% Interdiction de coupures de mots en fin de page
\brokenpenalty = 10000 

% Table des matières ---------------------------------------------------------------

\renewcommand\tableofcontents{{%
    \if@twocolumn
      \@restonecoltrue\onecolumn
    \else
      \@restonecolfalse
    \fi
    \parskip \z@%
    \iftreatise \setcounter{page}{9}
    \else          \setcounter{page}{5}
    \fi
    \chapter*{\contentsname}%
    \@starttoc{toc}%
    \if@restonecol\twocolumn\fi
    }}

\setcounter{tocdepth}{3}   % On place normalement dans la table des matières jusqu'aux
                           % intertitres de niveau 3 (réduire si nécessaire pour avoir
                           % un nb total de pages multiple de 16 ou 32

\renewcommand{\@pnumwidth}{6.35mm} % largeur de la colonne de droite pour les no pages
\renewcommand{\@tocrmarg}{10mm} % largeur marge de droite pour les titres sur plusieurs
                                % lignes

\def\TextInToc#1#2#3{%
% #1 Niveau de l'unité de sectionnement
% #2 indentation gauche (la marge de droite est à \@tocrmarg du bord droit du texte
% #3 texte complet à insérer dans la table des matières, par exemple nom d'auteur
  \ifnum #1>\c@tocdepth \else %
    \vskip \z@ \@plus.2\p@ %
    {\leftskip #2\relax \rightskip \@tocrmarg \parfillskip -\rightskip
     \parindent #2\relax\@afterindenttrue
     \interlinepenalty\@M
     \leavevmode
     {#3}\nobreak \hfill \nobreak \hb@xt@\@pnumwidth{\empty}%
    \par}
  \fi
}

\newcommand{\AddToContents}[3]{%
% Éventuellement, pour l'utilisateur, pour des cas spéciaux
% #1 Niveau de l'unité de sectionnement
% #2 indentation gauche (la marge de droite est à \@tocrmarg du bord droit du texte
% #3 texte complet à insérer dans la table des matières, par exemple nom d'auteur
\addtocontents{toc}{\protect\TextInToc{#1}{#2}{#3}}}

% espacement entre les points de suite en mu (math units) pardéfaut dans book : 4.5
\renewcommand{\@dotsep}{4}  % correspond expérimentalement à un espace de 2mm entre 2 points

% Lignes de sections dans la table des matiéres

\def\@dottedfiletoc{% D'après \@dottedtoclinedans latex.ltx
    % Nécessaire pour \l@part, \l@chap pour contrôler
    % l'espacement entre les points de \@dotsep (\dotfill est un peu différent)
    % points de suite (espacement de \@dotsep mu)
                         %\fontsize{8pt}{12pt}\selectfont très fin
   {\FonteTexte\mdseries\leaders\hbox{$\m@th
   \mkern \@dotsep mu\hbox{.}\mkern \@dotsep
   mu$}\hfill}}

\renewcommand*\l@part[2]{% Adaptation de book.cls
% #1 texte complet de \contentsline SANS \numberline{...} de la forme :
%    PREMIÈRE PARTIE. Titre de la partie
% #2 numéro de page
  \ifnum \c@tocdepth >-2\relax
    \addpenalty{-\@highpenalty}         % Autorise un saut de page
    \addvspace{11pt plus 1pt minus 1pt} % Espace vertical au dessus
     {\leftskip 0pt\relax               % Marge gauche
      \rightskip \@tocrmarg        % Marge droite (si plusieurs lignes, très déconseillé
      \parfillskip -\rightskip     % justification jusqu'à marge droite
     \parindent 0pt\relax\@afterindenttrue 
     \interlinepenalty\@M          % pénalité moyenne de cupure de page entre les lignes
     \leavevmode                   % passe en mode horizontal
     \textbf{#1}                   % Titre en gras
     \nobreak\@dottedfiletoc       % points de suite en maigre (espacement de \@dotsep mu)
     \nobreak\hb@xt@\@pnumwidth{\hss #2} % numéro de page cadré à droite ds boîte de \@pnumwidth
     \par}
     \penalty\@highpenalty         % Interdit un saut de page
  \fi
}

\newlength{\@toclmargchap}
\setlength{\@toclmargchap}{0pt}

\renewcommand*\l@chapter[2]{% Adaptation de book.cls
% #1 texte complet de \contentsline SANS \numberline{...} de la forme :
%    Chapitre 1. Titre de chapitre
% #2 numéro de page
  \ifnum \c@tocdepth >\m@ne
    \addpenalty{-\@highpenalty}% Autorise un saut de page
    \vskip 11pt \@plus1\p@ %N au lieu de 1em
     {\leftskip \@toclmargchap \relax  % Marge gauche
      \rightskip \@tocrmarg        % Marge droite (si plusieurs lignes, très déconseillé
      \parfillskip -\rightskip     % justification jusqu'à marge droite
     \parindent \@toclmargchap\relax\@afterindenttrue 
     \interlinepenalty\@M          % pénalité moyenne de cupure de page entre les lignes
     \leavevmode                   % passe en mode horizontal
     \textbf{#1}                   % Titre en gras
     \nobreak\@dottedfiletoc       % points de suite en maigre
     \nobreak\hb@xt@\@pnumwidth{\hss #2} % numéro de page cadré à droite ds boîte de \@pnumwidth
     \par}
     \penalty\@highpenalty         % Interdit un saut de page
  \fi
}

% Lignes de sections dans la table des matiéres

\def\@dottedtocline#1#2#3#4#5{% Adaptation de latex.txt pour appeler \@dottedfiletoc,
                              % une adaptation du tracé des points de suite.
  \ifnum #1>\c@tocdepth \else
    \vskip \z@ \@plus.2\p@
    {\leftskip #2\relax \rightskip \@tocrmarg \parfillskip -\rightskip
     \parindent #2\relax\@afterindenttrue
     \interlinepenalty\@M
     \leavevmode
     \@tempdima #3\relax
     \advance\leftskip \@tempdima \null\nobreak\hskip -\leftskip
     {#4}%
     \nobreak\@dottedfiletoc % points de suite en maigre
     \nobreak                %
     \hb@xt@\@pnumwidth{\hfil\normalfont \normalcolor #5}%
     \par}%
  \fi}

\renewcommand{\l@section}{\@dottedtocline{1}{5mm}{6mm}}
\renewcommand{\l@subsection}{\@dottedtocline{2}{10mm}{8.5mm}}
\renewcommand{\l@subsubsection}{\@dottedtocline{3}{15mm}{11mm}}

% Lignes d'annexes dans la table des matiéres

\def\l@appendix{\@dottedtocline{1}{5mm}{4mm}}
\def\l@appendixi{\@dottedtocline{1}{10mm}{7mm}}
\def\l@appendixii{\@dottedtocline{2}{15mm}{9.5mm}}

% Notes de bas de page -----------------------------------------------------------------

% corps 9pt/11pt, numérotation décimale
% séparées du texte par 1 ligne blanche
% filet maigre de 2.5 cm

% Espace au dessus du filet
\setlength{\skip\footins}{7mm} % mesuré sur Exemple fourni

\global\newcount\@prevpage     % pour détecter les footnotes à cheval sur 2 pages

\long\def\@footnotetext#1{%
\protect\global\@prevpage\c@page %NEW
    \insert\footins{%
    \reset@font\footnotesize
    \interlinepenalty\interfootnotelinepenalty
    \splittopskip\footnotesep
    \splitmaxdepth \dp\strutbox \floatingpenalty \@MM
    \hsize\columnwidth \@parboxrestore
    \protected@edef\@currentlabel{%
       \csname p@footnote\endcsname\@thefnmark
    }%
    \color@begingroup
      \@makefntext{%
        \rule\z@\footnotesep\ignorespaces#1\@finalstrut\strutbox}%
    \color@endgroup}}%

\long\def\@sfootnotetext#1{%
   % Pour faire une footnote sans numéro
   % = \@footnotetext de latex.ltx avec appel à \@makesfntext
    \protect\global\@prevpage\c@page %NEW
    \insert\footins{%
    \interlinepenalty\interfootnotelinepenalty
    \splittopskip\footnotesep
    \splitmaxdepth \dp\strutbox \floatingpenalty \@MM
    \hsize\columnwidth \@parboxrestore
    \protected@edef\@currentlabel{%
       \csname p@footnote\endcsname\@thefnmark
    }%
    \color@begingroup
      \@makesfntext{%NEW
        \rule\z@\footnotesep\ignorespaces#1\@finalstrut\strutbox}%
    \color@endgroup}}%

% Je n'ai pas réussi à automatiser la définition de \footnoterule lorsque c'est une continuation
\gdef\largefootnoterule{\par\vspace*{-3.3pt}\noindent\rule{12cm}{0.3pt}\vspace*{3pt}}
\gdef\smallfootnoterule{\par\vspace*{-3.3pt}\noindent\rule{2.5cm}{0.3pt}\vspace*{3pt}}
\newcommand{\LargeFnRule}{\global\let\footnoterule\largefootnoterule}
\newcommand{\SmallFnRule}{\global\let\footnoterule\smallfootnoterule}
\SmallFnRule

% Espace vertical avant chaque note (strut)
% Si \footnotesep est < \baselineskip (11pt), cela n'a pas d'effet
\setlength{\footnotesep}{9pt} % expérimental

% Appel de note : 8pts, élevé de 3 pts au dessus de la ligne
% Il vaut mieux utiliser des définitions relatives à la police courante
% que 3pt et 8 pt

\def\@makefnmark{\raisebox{0.55ex}{\FonteAppelNote \@thefnmark}}

\def\@Myfnmark{% Redéfinit la marque de footnote comme un nombre ordinaire
      % french.sty redéfinit bêtement \@makefntext => donc j'utilise une autre macro.
      \mbox{\arabic{footnote}.~}%
}

\renewcommand{\@makefntext}[1]{%
   % Place une footnote en bas de page avec le numéro \thefootnote
   % \@prevpage sert à détecter les notes de bas de page à cheval sur 2 pages
   \FonteNoteBasPage %NEW
   \noindent\@Myfnmark#1 %
}%

\def\@makesfntext#1{
   % Place une footnote sans numéro, pour les noms d'auteurs par exemple
   \FonteNoteBasPage %NEW
   \noindent #1 %
}%

\def\@thefnmark{\normalsize\arabic{footnote}}

% Objets flottants: figures, tableaux et programmes -------------------------------------

% Correction du bug du style french qui redéfinit le corps des légendes par \em
% et non \itshape (gênant si l'utilisateur redéfinit \em).

% Le plus fiable est de redéfinir \@caption de latex.ltx qui est l'unique appel
% à \@makecaption (En espérant que \@caption ne sera pas redéfini dans latex).

\def\Center{% redéfinition de l'environnement center sans espaces verticaux supplémentaires
  \let\@trivlistORI\@trivlist%
  \def\@trivlist{%
  \itemsep=0pt\parsep=0pt\topsep=0pt\partopsep=0pt\parskip=0pt%
  \@trivlistORI}%
  \trivlist \centering\item\relax}
\def\endCenter{%
  \endtrivlist}

\newlength{\rulewidth}\setlength{\rulewidth}{0.25pt}
\setlength{\arrayrulewidth}{\rulewidth}
\setlength{\fboxrule}{\rulewidth}
\setlength{\fboxsep}{0pt}

\long\def\@caption#1[#2]#3{{%
   \parskip \z@\par
% @@ revoir la mise des légendes dans la ToC, car il faut un argument facultatif de légende
%   courte (sans coupures de lignes) pour traiter cela correctement.
%   \addcontentsline{\csname ext@#1\endcsname}%
%                   {#1}{\protect\numberline{\csname the#1\endcsname}{\ignorespaces #2}}%
   \begingroup%
     \@parboxrestore%
     \normalsize%
     \@Mymakecaption{\csname fnum@#1\endcsname}{\ignorespaces #3}%
     % L'espace sous la légende est obtenu par \textDfloatsep ou \textfloatsep (cf ci-dessous)
   \endgroup}}

\newdimen\@@vspacapt \@@vspacapt=14.2pt
%\settoheight{\@@vspacapt}{\FonteTypeLegende T} = 6.2pt + 8 = 14.2pt

\newdimen\captionwidth
\newif\if@bigcaption\@bigcaptionfalse

\def\@Mymakecaption#1#2{%
    % #1 est e.g. Figure 1 ou Tableau 1, #2 est la légende de la figure
   \if@bigcaption
   \parbox[t]{12cm}{%
      \FonteTypeLegende\rule{0pt}{14.2pt}% 8 pts au dessus
       #1.\hspace{3pt}\FonteTitreLegende #2}
   \else
   \parbox[t]{12cm}{%
      \centering\FonteTypeLegende\rule{0pt}{14.2pt}% 8 pts au dessus
       #1.\hspace{3pt}\FonteTitreLegende #2}
   \fi
   }%

% Tableaux et Figures ------------------------------

% On définit l'environnement Table pour fixer les dimensions requises
% et pour placer la légende où il faut.
\setlength{\arrayrulewidth}{0.25pt}
\setlength{\tabcolsep}{2pt}
\def\@@caption{}

\newenvironment{Figure}[2][]{%
% #1 : placement de la figure du type htb... comme pour figure ou table
% #2 : la légende
\settowidth{\captionwidth}{\FonteTitreLegende #2}
%\renewcommand{\@@caption}{#2\\(~\the\captionwidth~)}%
\renewcommand{\@@caption}{#2}%
\ifdim\captionwidth > 345mm \@bigcaptiontrue \else\@bigcaptionfalse\fi
\@float{figure}[#1]%
\begin{Center}
}{\caption{\@@caption}\end{Center}\end@float}

% Pour les utilisateurs peu familiers avec epsfig, nous fournissons ceci :
\newcommand{\EpsFigure}[4][]{%
% #1 : placement de la figure du type htb... comme pour figure ou table
% #2 : la légende
% #3 : la hauteur de la figure (par défaut, les dimensions d'origine)
% #4 : le fichier qui contient la figure en format .eps
\begin{Figure}[#1]{#2} %
\ifthenelse{\equal{#3}{}}{\epsfig{file=#4}}{\epsfig{height=#3,file=#4}} %
\end{Figure} %
}

\newenvironment{Table}[2][]{%
% #1 : placement de la table du type htb... comme pour figure ou table
% #2 : la légende
\settowidth{\captionwidth}{\FonteTitreLegende #2}
\renewcommand{\@@caption}{#2}%
\ifdim\captionwidth > 345mm \@bigcaptiontrue \else\@bigcaptionfalse\fi
\@float{table}[#1]%
\renewcommand{\arraystretch}{1} % modifie l'espace entre les lignes de texte
\begin{Center}
\FonteNoteBasPage %
}{\caption{\@@caption}\end{Center}\end@float}

% Flottants -------------------------------------------------------------------------

\def\topfigrule{}               % cmde de tracé 1 règle float haut et texte
\def\botfigrule{}               % cmde de tracé 1 règle float bas et texte

% Ajustement des paramètres des figures flottantes
% L'utilisateur PEUT, voire DOIT AJUSTER ces valeurs s'il rencontre
% des problèmes de placements ! cf Latex companion

\setcounter{topnumber}{2}              % déf. 2    % nb max de fig. en haut
\setcounter{bottomnumber}{1}           % déf. 1    % nb max de fig. en bas
\setcounter{totalnumber}{3}            % déf. 3    % nb max de fig. par page
\renewcommand{\topfraction}{1}         % déf. 0.7  % fraction max de remplissage de la
                                                   % page par la (les) figure(s) en haut
\renewcommand{\bottomfraction}{1}      % déf. 0.3  % fraction max de remplissage de la
                                                   % page par la figure en bas
                                                   % en baspar fig
\renewcommand{\textfraction}{0}        % déf. 0.2  % fract. min remplissage texte / page
\renewcommand{\floatpagefraction}{0.8} % déf. 0.5  % remplissage floatpage

% Dimensions des espaces verticaux élastiques autour des figures :
% jouent un rôle important pour le placement, on a intérêt à augmenter l'élasticité

%\setlength{\floatsep}{22pt plus 5pt minus 2pt}  
\setlength{\floatsep}{22pt}  
       % espace supplémentaire entre 2 figures
       % par défaut: \setlength{\floatsep}{12pt plus 2pt minus 2pt}

%\setlength{\textfloatsep}{22pt plus 2pt minus 4pt}
\setlength{\textfloatsep}{22pt}
       % espace entre texte et figure pour des figures en haut ou en bas
       % par défaut: \setlength{\textfloatsep}{20pt plus 2pt minus 4pt}

\newskip\intextUsep
\setlength{\intextUsep}{\parskip} 
       % espace supplémentaire au dessus d'une figure [h], entre figure et texte.
       % par défaut: \setlength{\intextsep}{\floatsep}

\newskip\intextDsep
\setlength{\intextDsep}{22pt plus 5pt minus 2pt} 
       % espace supplémentaire au dessous d'une
       % figure [h], entre figure et texte.
       % par défaut: \setlength{\intextsep}{\floatsep}

% redéfinition de Latex.ltx pour avoir un espace vertical au dessus et
% au dessous des figures flottantes au milieu différent (par défaut c'est le même au
% dessus et au dessous : \intextsep
\def \@addtocurcol {% adaptation de latex.ltx
   \@insertfalse
   \@setfloattypecounts
   \ifnum \@fpstype=8
   \else
     \ifnum \@fpstype=24
     \else
       \@flsettextmin
       \advance \@textmin \@textfloatsheight
       \@reqcolroom \@pageht
       \ifdim \@textmin>\@reqcolroom
         \@reqcolroom \@textmin
       \fi
       \advance \@reqcolroom \ht\@currbox
       \ifdim \@colroom>\@reqcolroom
         \@flsetnum \@colnum
         \ifnum \@colnum>\z@
           \@bitor\@currtype\@deferlist
           \if@test
           \else
             \@bitor\@currtype\@botlist
             \if@test
               \@addtobot
             \else
               \ifodd \count\@currbox
                 \advance \@reqcolroom \intextUsep          %%% 1er
                 \ifdim \@colroom>\@reqcolroom
                   \global \advance \@colnum \m@ne
                   \global \advance \@textfloatsheight \ht\@currbox
                   \global \advance \@textfloatsheight \intextUsep     %%% 1er 
                   \global \advance \@textfloatsheight \intextDsep     %%% + 2eme
                   \@cons \@midlist \@currbox
                   \if@nobreak
                     \nobreak
                     \@nobreakfalse
                     \everypar{}%
                   \else
                     \addpenalty \interlinepenalty
                   \fi
                   \vskip \intextUsep      %%% 1er
                   \box\@currbox
                   \penalty\interlinepenalty
                   \vskip\intextDsep       %%% 2eme
                   \ifnum\outputpenalty <-\@Mii \vskip -\parskip\fi
                   \outputpenalty \z@
                   \@inserttrue
                 \fi
               \fi
               \if@insert
               \else
                 \@addtotoporbot
               \fi
             \fi
           \fi
         \fi
       \fi
     \fi
   \fi
   \if@insert
   \else
     \@resethfps
     \@cons\@deferlist\@currbox
   \fi
}
% Formules -------------------------------------------------------------------------------

% En mode standard, les modes equation, eqnarray et eqnarray*
% fournissent bien une ligne blanche au dessus et au dessous et la numerotation
% collée sur la marge de droite.
% La seule adadaptation est la numérotation [n] au lieu de (n).

\renewcommand{\@eqnnum}{\reset@font\rmfamily [\theequation]}

% Définition des espaces sans l'option fleqn :

\setlength{\mathindent}{10mm}           % indentation des formules
\setlength{\arraycolsep}{5pt}           % moitié de l'espace horizontal entre les colonnes
                                        % d'un environnement array
\setlength{\jot}{0.5\baselineskip}      % extra espace entre les lignes de eqnarray/eqnarray*
\setlength{\topsep}{11pt plus 3pt minus 9pt} % espace supplém. au dessus/dessous des formules.

% Les valeurs suivantes ne sont pas utilisées avec l'option fleqn
%\setlength{\abovedisplayskip}{0mm}      % espace supplémentaire au dessus des longues formules
%\setlength{\abovedisplayshortskip}{0mm} % idem, formules courtes
%\setlength{\belowdisplayskip}{0mm}      % espace supplémentaire au dessous des longues formules
%\setlength{\belowdisplayshortskip}{0mm} % idem, formules courtes

% Utilitaires facultatifs ---------------------------------------------------------------

% Utiliser la commande \eqncont pour terminer une ligne de formule
% sur plusieurs lignes : supprime la numérotation + espace vertical plus serré.
\newcommand{\EqnCont}{\nonumber \\[0\baselineskip]}

%%% Retardement du pagebreak

%  Pour retarder le passage automatique à la page suivante de 1 ou 2 lignes,
%  utiliser le couple de commandes \delaynewpage{nb lignes} ... \forcenewpage
%  Attention! cette technique est du bricolage, mais cependant indispensable
%  pour une présentation optimale d'un document. (indiquée par Lamport).
%  IL NE FAUT L'EMPLOYER QUE TOUT À LA FIN, JUSTE AVANT L'IMPRESSION FINALE.
%  En principe ne pas retarder de plus de 2 lignes !

\newcommand{\Delaynewpage}[1]{\enlargethispage{#1\baselineskip}}
   % À utiliser sur la page à allonger, à un endroit bien visible
   % avant le premier paragraphe ou la figure à faire tenir sur la page.
   % Provoque une descente des notes de bas de page s'il y en a !

\newcommand{\Forcenewpage}{\pagebreak\noindent}
   % À utiliser à l'endroit précis où la coupure doit être forcée.
   % on peut placer cette commande en début de ligne, avec suite du texte
   % sur la ligne suivante pour faciliter son repérage

% Listes ----------------------------------------------------------------------------

% Espaces communs aux listes numérotées (enumerate) ou non (itemize)
% Attention !!! Ces espacements ne peuvent être efectifs que si Babel ne les redéfinit
% pas dans \@trivlist, comme il tente de l'imposer. Il faut ABSOLUMENT
% exécuter \bbl@nonfrenchlistspacing pour revenir à l'espacement Hermes,
% en particulier après \selectlanguage qui redéfinit l'espacement de babel.

\def\@hermesitemizespacing{%
% espacements horizontaux
   \leftmargin=0pt
   \labelsep=1mm
   \ifnum \@itemdepth < 2  % niveau 0 ou 1
      \setlength{\itemindent}{\itemindenti}
      \listparindent=5mm
   \else                   % niveau 2
      \itemindent=\itemindentii
      \listparindent=10mm
   \fi
   \labelwidth=0pt
   \rightmargin=0pt
  % espacements verticaux des listes
    \ifnum \@itemdepth = 0 %
         \itemsep=0pt   %  niveau 0
         \topsep=3pt 
    \else\ifnum \@itemdepth = 1 %
         \parskip=0pt
         \itemsep=3pt %  niveau 1
         \topsep=3pt
    \else
         \topsep=0pt
         \parskip=0pt
         \itemsep=0pt %  niveau 2
    \fi\fi
    \partopsep=0pt
    \parsep=0pt
% pénalités
   \sloppy%
   \widowpenalty10000%ligne isolée en début de page
   \clubpenalty10000%ligne isolée à la fin d'une page
}

\newcommand{\displaylistdimen}{% 
  % Utilitaire de mise au point en cas de réapparition du problème de l'espacement de
  % babel : cf hermesitemize et hermesenumerate
     itemdepth: \the\@itemdepth ;\\
     topsep: \the\topsep ;\\
     itemsep: \the\itemsep ;\\
     parskip: \the\parskip ;\\
     partopsep: \the\partopsep ;\\
     parsep: \the\parsep.
}

% Listes non numérotées (itemize) -----------------------------------------------------------------------

\newlength{\labelwidthi}      % largeur de la boîte contenant le retrait de 5mm + tiret moyen
\settowidth{\labelwidthi}{--} % environ 2mm, tiret moyen
\addtolength{\labelwidthi}{5mm} % retrait + largeur du tiret moyen, environ 7mm

\newlength{\itemindenti}      % décalage du texte de 1ère ligne depuis marge gauche 
\setlength{\itemindenti}{\labelwidthi} % = labelwidthi
\addtolength{\itemindenti}{1mm}          %   + labelsep

\newlength{\labelwidthii}       % largeur de la boîte contenant le retrait de 10mm + tiret court
\settowidth{\labelwidthii}{-}   % environ 1mm, tiret court
\addtolength{\labelwidthii}{10mm} % retrait + largeur du tiret court, environ 11mm

\newlength{\itemindentii}      % décalage du texte de 1ère ligne depuis marge gauche 
\setlength{\itemindentii}{\labelwidthii} % = labelwidthii
\addtolength{\itemindentii}{1mm}           %   + labelsep

\newenvironment{hermesitemize}{% liste à deux niveaux seulement
 \ifnum \@itemdepth > 1 \@toodeep % ici, \@itemdepth = niveau -1 
 \else%
    \renewcommand{\labelitemi}{{\makebox[\labelwidthi][r]{--}}}%
    \renewcommand{\labelitemii}{{\makebox[\labelwidthii][r]{-}}}%
    \advance\@itemdepth \@ne%
    \edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
    \begin{list}%
       {\csname\@itemitem\endcsname}%
       {\def\makelabel##1{\llap{##1}}% ##1 signifie arg 1 de list et non de hermesitemize
        % définir les espaces ici et nulle part ailleursn
        \@hermesitemizespacing %
       }%\item \displaylistdimen
 \fi}{\end{list}}%

% Listes numérotées ---------------------------------------------------------------------------------------

\newlength{\enumwidthi}      % largeur de la boîte contenant le retrait de 5mm + n.
\settowidth{\enumwidthi}{5)} % environ 2mm, chiffre le plus large
\addtolength{\enumwidthi}{5mm} % retrait + largeur du tiret moyen, environ 7mm

\newlength{\enumindenti}      % décalage du texte de 1ère ligne depuis marge gauche 
\setlength{\enumindenti}{\enumwidthi} % = enumwidthi
\addtolength{\enumindenti}{1mm}          %   + enumsep

\newlength{\enumwidthii}         % largeur de la boîte contenant le retrait de 10mm + tiret court
\settowidth{\enumwidthii}{b)}   % environ 1mm, tiret court
\addtolength{\enumwidthii}{10mm} % retrait + largeur du tiret court, environ 11mm

\newlength{\enumindentii}      % décalage du texte de 1ère ligne depuis marge gauche 
\setlength{\enumindentii}{\enumwidthii} % = enumwidthii
\addtolength{\enumindentii}{1mm}           %   + enumsep

\newenvironment{hermesenumerate}{% enumération à un niveau
 \ifnum \@itemdepth >1 \@toodeep
 \else%
      \advance\@enumdepth \@ne%
      \advance\@itemdepth \@ne%
      \edef\@enumctr{enum\romannumeral\the\@enumdepth}%
      \begin{list}%
         {\csname label\@enumctr\endcsname}%
         {%
          \usecounter{\@enumctr}%
          \def\makelabel##1{\llap{##1}}% ##1 signifie arg 1 de list et non de hermesenumerate
          % définir les espaces horizontaux ici
          \@hermesitemizespacing%
          \ifnum \@itemdepth =1
             \setlength{\itemindent}{\enumindenti}
          \else
             \setlength{\itemindent}{\enumindentii}
          \fi
          \renewcommand{\labelenumi}{{\makebox[\labelwidthi][r]{\arabic{enumi})}}}
          \renewcommand{\labelenumii}{{\makebox[\labelwidthii][r]{\alph{enumii})}}}%
         }
   \fi}{\end{list}}


% Quotations, verses... redéfinition des espaces verticaux à 6 pts

\renewenvironment{quote}{%
\begin{list}{%
  \itemsep=0pt\parsep=0pt\topsep=0pt\partopsep=3pt\parskip=3pt%
  \leftmargin=1cm\labelsep=0pt\listparindent=0pt\itemindent=0pt\labelwidth=0pt\rightmargin=0pt%
  }\item\relax}{\end{list}}

% Remarques, notes, exemples, Définitions, Théorème... ----------------------------------------------------

\newcommand{\GenericRemark}[2]{{%
  % #1 : type de remarque, par exemple Note, Exemple, ...
  % #2 : contenu de la remarque
    \par\parindent 0pt\textsc{#1}.--\ \nopagebreak[4] {\parskip 5pt\relax#2}}}

\newcommand{\Remark}[1]{\GenericRemark{\remarkname}{#1}}
\newcommand{\Note}[1]{\GenericRemark{\notename}{#1}}
\newcommand{\Example}[1]{\GenericRemark{\examplename}{#1}}

% Cas des remarques sur plusieurs paragraphes
\newenvironment{Remarks}{\setlength{\parskip}{5pt}%
    \par\noindent\textsc{\remarkname}.--\ \nopagebreak[4]}{\setlength{\parskip}{\@normalparskip}}

\newenvironment{Notes}{\setlength{\parskip}{5pt}%
    \par\noindent\textsc{\notename}.--\ \nopagebreak[4]}{\setlength{\parskip}{\@normalparskip}}

\newenvironment{Examples}{\setlength{\parskip}{5pt}%
    \par\noindent\textsc{\examplename}.--\ \nopagebreak[4] }{\setlength{\parskip}{\@normalparskip}}

\newcommand{\GenericDefinition}[2]{{%
  % #1 : type de definition, par exemple Définition, Théorème, ...
  % #2 : contenu de la définition
    \par\parindent 0pt\textsc{#1}.--\ \nopagebreak[4] {\parskip 5pt\relax\textit{#2}}}}

\newcommand{\Definition}[1]{\GenericDefinition{\definitionname}{#1}}
\newcommand{\Theorem}[1]{\GenericDefinition{\theoremname}{#1}}
\newcommand{\Lemma}[1]{\GenericDefinition{\lemmaname}{#1}}

% Bibliographie --------------------------------------------------------------------------

% J'ai séparé la définition de thebibliography en 2 parties
% thebibliography ajoute le nom de section puis appelle la liste BiblioHermes

\iftreatise
\else \bibliographystyle{ouvrage-hermes}
\fi

\renewenvironment{thebibliography}[1]{%
   \@endofappendices
   \iftreatise \section{\bibname}
   \else \chapter{\bibname} \fi
   \begin{BiblioHermes}%
  }{%
   \end{BiblioHermes}%
}

% BiblioHermes affiche une liste de references bibliographiques.
% Les entrées sont normalement préparées par le style ouvrage-hermes.bst, mais
% peuvent aussi s'écrire à la main :
% \bibitem[ref]{cite} qui écrit comme d'habitude
% la clef de citation dans le fichier .aux pour bibtex
% \bibi[ref] qui affiche seulement l'entree, sans écrire sur .aux
% (utilité très spéciale, pour placer plusieurs biblios dans un document).

\def\bibi[#1]{\item[\@biblabel{#1}\hfill]} % @ special

\newenvironment{BiblioHermes}{%begin
   \begin{list}{}{% Définition des paramètres de liste pour la bibliographie Hermes
         \usecounter{enumiv}%
         \let\p@enumiv\@empty %
         \renewcommand\theenumiv{\arabic{enumiv}}%
         \renewcommand\newblock{\hskip .11em \@plus.33em \@minus.07em}%
         \FonteBibliographie%
         %% dimensions horizontales
         \setlength{\leftmargin}{5mm}%
         \setlength{\itemindent}{-3mm}%
         \setlength{\labelsep}{2mm}%
         \setlength{\labelwidth}{0mm}%
         %% dimensions verticales
          \setlength{\topsep}{0pt}%
          \setlength{\parskip}{11pt plus 2.2pt minus 1.1pt}%
          \setlength{\itemsep}{5pt plus 1pt minus .5pt}%
          \setlength{\partopsep}{0pt}%
          \setlength{\parsep}{0pt}%
         \sloppy\clubpenalty10000\widowpenalty10000%
         \sfcode`\.=\@m %
   }}{
      \def\@noitemerr{\@latex@warning{Empty `thebibliography' environment}}%
      \FonteTexte%
      \end{list}}

\renewcommand{\numberline}[1]{#1}

% Index --------------------------------------------------------------------------------

% L'index doit être produit avec makeindex et le fichier de style ouvrage-hermes.ist
% (cf le manuel de makeindex sur votre machine pour savoir comment spécifier ce fichier
% de style), par exemple :
% makeindex -l -s ouvrage-hermes.ist mon_livre.idx % pour un index unique
% makeindex -l -s ouvrage-hermes.ist mon_livre_chap1.idx ... mon_livre_chap2.idx 
            % pour des index par chapitres

% Les macros suivantes sont pour la plupart des redéfinitions de book.cls

\renewenvironment{theindex}{%
  \if@openright\cleardoublepage\else\clearpage\fi %
    \columnseprule \z@ % 0.25pt
    \columnsep 1cm %
    \begin{multicols}{2}[\chapter{\indexname} %chapitre non numéroté puisque index dans backmatter
           \enlargethispage{0pt}] %correction du décalage inexpliqué sur 1ere page
    \begin{flushleft}%
    \FonteNoteBasPage %
    \parindent\z@ %
    \parskip\z@ \@plus .3\p@\relax %
    \let\item\@idxitem%
  }{\end{flushleft}\end{multicols}}

%\par\hangindent x mm fixe la marge gauche des lignes de continuation
%le décalage éventuel de chaque item est obtenu par un hspace

\renewcommand\@idxitem{\par\hangindent 5mm} 
   % entrée niveau 1, sans alinéa,
   % avec retrait de 5mm sur ligne de continuation

\renewcommand\subitem{\par\hangindent 10mm\hspace*{5mm}}
   % entrée de niveau 2, avec alinéa de 5mm,
   % et retrait de 10mm sur ligne de continuation

\renewcommand\subsubitem{\par\hangindent 15mm \hspace*{10mm}} 
   % entrée de niveau 3, avec alinéa de 10mm,
   % et retrait de 15mm sur ligne de continuation

\renewcommand\indexspace{}
   % espace avant lettre de séparation

\def\mkidxletter#1{%
   \ifthenelse{\equal{#1}{A}}{}{\par \vskip 11\p@ \@plus3\p@ \@minus2\p@\relax}
   \textbf{#1}\nopagebreak\par \vskip 5\p@\relax}
   % place la lettre de séparation en gras, fer à gauche.
   % suivie d'un espace vertical de 5pts (spécifique de ouvrage-hermes)

% Fiche pour le service de fabrication -------------------------------------------------
\def\@@soustitre#1{{\parskip3pt\par\noindent\textbf{#1}}\par}

\newcommand{\Publisher}[3]{
   % Pour construire la page d'annexe pour le service de fabrication
   %        {#1} = numero de téléphone des auteurs
   %        {#2} = fax des auteurs
   %        {#3} = Email des auteurs
   \renewcommand{\VAR@Phone}{#1}
   \renewcommand{\VAR@Fax}{#2}
   \renewcommand{\VAR@Email}{#3}
   \AtEndDocument { { \setlength{\parskip}{6pt}
       \clearpage
       \raggedbottom
       \vspace*{10pt}
       \begin{center}
           \iffrenchlang \textbf{\Large Fiche pour le service de fabrication}
           \else          \textbf{\Large Note to Publisher}
           \fi
       \end{center}
       \vspace*{5pt}
       \section*{\iffrenchlang Auteurs~: \else Authors: \fi}
         \begin{quote} \@VARAuteurs \end{quote}

      \section*{\iffrenchlang Titre du livre~: \else Book title: \fi}
         \begin{quote}\parbox{6cm}{\raggedright \@VARTitreComplet } \end{quote}

      \section*{\iffrenchlang Titre abrégé~: \else Short title: \fi}
         \begin{quote}\textit{\VAR@TitreAbrege}\end{quote}

      \section*{\iffrenchlang Date de cette version~: \else Date of this version: \fi}
         \begin{quote}\textit{\today}\end{quote}

      \section*{\iffrenchlang Contact~: \else Contact: \fi}
        \vspace{6pt}
        \begin{itemize}
             \item \iffrenchlang téléphone~: \else phone: \fi \VAR@Phone
             \item \iffrenchlang télécopie~: \else fax: \fi   \VAR@Fax 
             \item \iffrenchlang courriel~: \else email: \fi \VAR@Email
        \end{itemize}

      \section*{\iffrenchlang Logiciel pour la composition~: \else Typesetting system: \fi}
       \vspace{6pt}
      \begin{itemize}
        \item \LaTeX, \iffrenchlang avec la classe\else with the class \fi
              \texttt{\PackageName.cls}
        \item version~\PackageVersion, \filedate
        \item \iffrenchlang 
                 traité (option treatise)~: 
                 \iftreatise Oui \textit{(chapitres avec différents auteurs)}
                 \else Non \textit{(par défaut, livre, mêmes auteurs de chapitres)} \fi
              \else
                 chapters with several authors (treatise option): 
                 \iftreatise Yes \else No \textit{(by default)} \fi
              \fi
        \item \iffrenchlang 
                livre en anglais (option english)~: Non \textit{(par défaut en français)} 
              \else 
                book in english (english option): Yes 
              \fi
        \item \iffrenchlang 
                 tracé des limites de page (option cropmarks)~: 
                 \ifcropmarks Oui \else Non \textit{(par défaut)} \fi
              \else
                 crop marks for trimming the pages (cropmarks option): 
                 \ifcropmarks Yes \else No \textit{(by default)} \fi
              \fi
        \item \iffrenchlang 
                   suppression des en-têtes de page (option empty)~: 
                   \ifemptyheaders      Oui \else   Non \textit{(par défaut)} \fi
               \else
                   no page headings (empty option): 
                   \ifemptyheaders      Yes \else   No \textit{(by default)} \fi
               \fi
        \item \iffrenchlang 
                   impression des pages blanches (option allpages)~: 
                   \ifallpages  Oui \else   Non \textit{(par défaut)} \fi
              \else
                   printing blank pages (allpages option): 
                   \ifallpages  Yes \else   No \textit{(by default)} \fi
              \fi
        \item \iffrenchlang 
                   césures actives~: voir la coupure du mot <<~cognitive~>> dans le fichier .log
              \else
                   active hyphenation: see hyphenation of the word ``cognitive'' in the .log file
              \fi
      \showhyphens{cognitive}
      \end{itemize}
      }}}% \publisher

%--- Internationalisation  Macros utilisées pour le paramétrage sur le langage
% Beaucoup de ces valeurs sont définies par babel, mais on ne sait jamais :
% les valeurs suivantes sont utilisées par ouvrage-hermes.cls ou ouvrage-hermes.bst

%--- Macros identiques en francais et english

\def\pagename{p.~}
\def\pagesname{p.~}
\def\etalter{\ et al.}
\def\Inname{}
\def\inname{}
\def\editornames{Eds.}
\def\editorname{Ed.}
\def\Volumename{Vol.~}
\def\volumename{vol.~}
\def\no{n°} % redéfinition de frenchb (préférence de Hermès)
\def\No{N°} % redéfinition de frenchb (préférence de Hermès)

\def\hermesup#1{\raisebox{0.55ex}{\protect\footnotesize #1\kern+.17em}}%
\def\hermesdown#1{\raisebox{-0.55ex}{\protect\footnotesize #1\kern+.17em}}%
\def\fsc{\scshape}%

%--- Macros différentes en francais et english
\def\@defand{%
   \iffrenchlang \def\and{ et } \else \def\and{ and } \fi
}
\@defand

\iffrenchlang
   \def\OrdinalInWords#1{\ifcase\value{#1} \or Première\or Deuxième\or Troisième
      \or Quatrième\or Cinquième\or Sixième\or Septième\or Huitième\or Neuvième
      \or Dixième\or Onzième\or Douzième\or Treizième\or Quatorzième\or Quinzième
      \or Seizième\or Dix-Septième\or Dix-Huitième\or Dix-Neuvième\or Vingtième\fi}%
\else
   \def\OrdinalInWords#1{\ifcase\value{#1} \or First\or Second\or Third
      \or Fourth\or Fifth\or Sixth\or Seventh\or Eighth\or Ninth
      \or Tenth\or Eleventh\or Twelfth\or Thirteenth\or Fourteenth\or Fifteenth
      \or Sixteenth\or Seventeenth\or Eighteenth\or Nineteenth\or Twentieth\fi}%
\fi

\def\@initialisation{%
  \iffrenchlang %
   \AutoSpaceBeforeFDP %
   \def\@partname{partie} %
   \def\chaptername{Chapitre}%
   \def\contentsname{Table des matières}%
   \def\listefigurename{Liste des figures}%
   \def\listetablename{Liste des tableaux}%
   \def\listeprogramname{Liste des programmes}%
   \def\figurename{Figure}%
   \def\figname{figure}%
   \def\tablename{Tableau}%
   \def\tabname{tableau}%
   \def\programname{Programme}%
   \def\progname{programme}%
   \def\HermesRefname{Bibliographie}%
   \def\abstractname{Résumé}%
   \def\remarkname{Remarque}%
   \def\notename{Note}%
   \def\examplename{Exemple}%
   \def\definitionname{Définition}%
   \def\theoremname{Théorème}%
   \def\lemmaname{Lemme}%
   \def\acknowledgementsname{Remerciements}%
   \def\appendixname{Annexe}%
   \def\appendicesname{Annexes}%
   \def\indexname{Index}%
   \def\bibname{Bibliographie}%
   \def\andname{, }%
   \def\editionname{édition}%
   \def\technicalreportname{Rapport}% valeur par défaut du champ type
   \def\numbername{\no}%{n°~} %\no est défini dans frenchb
   \def\Numbername{\No}%{N°~}
   \def\ofname{de}%
   \def\guilo{\og} % autre possibilité : \def\guilo{«~}
   \def\guilf{\fg} % autre possibilité : \def\guilf{~»}
   \def\Jan{janvier }%
   \def\Feb{février }%
   \def\Mar{mars }%
   \def\Apr{avril }%
   \def\May{mai }%
   \def\Jun{juin }%
   \def\Jul{juillet }%
   \def\Aug{août }%
   \def\Sep{septembre }%
   \def\Oct{octobre }%
   \def\Nov{novembre }%
   \def\Dec{décembre }%
\else %
   \gdef\@partname{part}%
   \def\chaptername{Chapter}%
   \def\contentsname{Contents}%
   \def\listefigurename{List of Figures}%
   \def\listetablename{List of Tables}%
   \def\listeprogramname{List of Programs}%
   \def\figurename{Figure}%
   \def\figname{figure}%
   \def\tablename{Table}%
   \def\tabname{table}%
   \def\programname{Program}%
   \def\progname{program}%
   \def\HermesRefname{References}%
   \def\abstractname{Summary}%
   \def\remarkname{Comment}%
   \def\notename{Note}%
   \def\examplename{Example}%
   \def\definitionname{Definition}%
   \def\theoremname{Theorem}%
   \def\lemmaname{Lemma}%
   \def\acknowledgementsname{Acknowledgements}%
   \def\appendixname{Appendix}%
   \def\appendicesname{Appendices}%
   \def\indexname{Index}%
   \def\bibname{Bibliography}%
   \def\andname{, }%
   \def\editionname{edition}%
   \def\technicalreportname{Report }% valeur par défaut du champ type
   \def\numbername{num.~}%
   \def\Numbername{Num.~}%
   \def\ofname{of}%
   \def\guilo{``}%
   \def\guilf{''}%
   \def\Jan{January }%
   \def\Feb{February }%
   \def\Mar{March }%
   \def\Apr{April }%
   \def\May{May }%
   \def\Jun{June }%
   \def\Jul{July }%
   \def\Aug{August }%
   \def\Sep{September }%
   \def\Oct{October }%
   \def\Nov{November }%
   \def\Dec{December }%
\fi%
}%\@initialisation

\def\@mkhermes{%
   \bbl@nonfrenchlistspacing%
   \gdef\itemize{\hermesitemize}\gdef\enditemize{\endhermesitemize}%
   \gdef\enumerate{\hermesenumerate}\gdef\endenumerate{\endhermesenumerate}%
   \let\up\hermesup%
   \let\down\hermesdown%
   \@initialisation % dépend du langage
}

\AtBeginDocument{
  \iffrenchlang %
      \selectlanguage{frenchb}%
  \else %
      \selectlanguage{english}%
      \let\@afterindentfalse\@afterindenttrue
      \@afterindenttrue
  \fi%
}

%EOF ouvrage-hermes.cls
