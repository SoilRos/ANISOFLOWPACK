% BibTeX bibliography style for Hermes publisher : 
% books, journals, proceedings, technical reports, thesis...
% R. Rousseau, rr@unice.fr, août 1999
% Remerciements : 
% Philippe.Queinnec@enseeiht.fr : MACRO {jan} {"\Jan{}"}...

% Utilitaires standards %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ENTRY
  { address
    author
    booktitle
    chapter
    edition
    editor
    howpublished
    institution
    journal
    key
    note
    number
    organization
    pages
    publisher
    school
    series
    title
    type
    volume
    month     % comme d'habitute, month et date concernent la date de publication
    year    
    date      % date correspond à la date d'une manifestation : la date de publication
              % doit être donnée aussi
  }
  {}
  { label extra.label sort.label }

INTEGERS { output.state before.all mid.sentence after.sentence after.block }

FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
}

STRINGS { s t }

FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.block =
        { % Hermes met des virgules partout => très laid
          % add.period$ write$
          ", " * write$
          newline$
          "\newblock " write$
        }
        { output.state before.all =
            'write$
            { % Hermes met des virgules partout => très laid
              % add.period$ " " * write$ 
              ", " * write$
            }
          if$
        }
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}

FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}

FUNCTION {output.guillemet}
{ duplicate$ empty$
    'pop$
    { "\guilo{}" swap$ * "\guilf{}" * output.nonnull }
  if$
}

FUNCTION {output.check} % S1 mess output.check
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}

FUNCTION {output.check.guillemet}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    { "\guilo{}" swap$ * "\guilf{}" * output.nonnull }
  if$
}

FUNCTION {output.bibitem}
{ newline$ 
  "\bibitem[" write$
  label write$
  "]{" write$
  cite$ write$
  "}" write$
  newline$
  ""
  before.all 'output.state :=
}

FUNCTION {fin.entry}
{ add.period$
  write$
  newline$
}

FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { after.block 'output.state := }
  if$
}

FUNCTION {new.sentence}
{ output.state after.block =
    'skip$
    { output.state before.all =
        'skip$
        { after.sentence 'output.state := }
      if$
    }
  if$
}

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}

FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}

FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}

FUNCTION {new.block.checka}
{ empty$
    'skip$
    'new.block
  if$
}

FUNCTION {new.block.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.block
  if$
}

FUNCTION {new.sentence.checka}
{ empty$
    'skip$
    'new.sentence
  if$
}

FUNCTION {new.sentence.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.sentence
  if$
}

FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}

FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "\textit{" swap$ * "}" * }
  if$
}

INTEGERS { nameptr namesleft numnames }

% Specifique du style biblio-hermes.bst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTION {format.names}
  % s        : STRING;              % chaîne argument à traiter
  % t        : STRING;              % chaîne résultat
  % nameptr  : INTEGER;             % indice du nom à traiter ( 1, 2, ...n)
  % nameleft : INTEGER;             % nombre de noms restant à traiter
  % numnames : INTEGER;             % nombre de noms à traiter au départ
{ 's :=                             % s := argument de la fonction
  #1 'nameptr :=                    % nameptr := 1
  s num.names$ 'numnames :=         % numnames := nombres de noms(s)
  numnames 'namesleft :=            % namesleft := numnames
  { namesleft #0 > }                % tantque ( namesleft > 0)
		  % format d'un groupe nominal :
		  % en petites capitales et  nom avant initiales des prénoms
        % rappel pour format.name$: 
        %   {v} particule {l} lastname {f} first name {j} junior
		  %   si on met 1 seule lettre, c'est l'initiale, avec 2 le nom entier.
   { 
      % t := groupe nominal [nameptr] formatté
      "\textsc{" s nameptr "{vv~}{ll}{ jj}{ f.}" format.name$ * "}" * 't :=
          nameptr #1 > % if indice du nom > 1      
          { % then ...
                namesleft #1 >   % il y a-t-il encore plus d'un nom à traiter?
	             { ", " * t * } % oui, séparateur ordinaire
	             { % cas du dernier nom 
                    % dernier séparateur s'il y a plus d'un nom en tout...
                    numnames #2 > % nombre total de noms > 2 ?
		              { "" * } % dernier séparateur si plus de 2 auteurs
                                   % On ne met pas de virgule: celle-ci sera eventuel.
                                   % placee dans \andname
		              'skip$        % dernier séparateur si 2 auteurs
	                 if$ % numnames #2 >
                    % si le dernier nom est "others" on écrit \etalter, sinon
                    % on intercale le mot \andname qui peut être redéfini par
                    % "et" ou "and"... dans le programme latex
	                   t "\textsc{others}" =
		                { "\etalter{}" * }
		                { "\andname{}" * t * }
	                 if$ % t "others" =
	             }
             if$ % namesleft #1 >
	       } % else
          't % indice du nom = 1, on l'écrit
      if$ % nameptr #1 >
      nameptr #1 + 'nameptr :=      % nameptr := nameptr + 1
      namesleft #1 - 'namesleft :=  % namesleft := namesleft - 1
   }
   while$
}

FUNCTION {format.authors}
% sortie du nom d'auteur, pour tous les types de documents
% On rajoute un terminateur : ","
{ author empty$
    { "" }
    { author format.names}
  if$
}

FUNCTION {format.editors}
{ editor empty$
    { "" }
    { editor format.names
      editor num.names$ #1 >
        { ", \editornames{}" * }
        { ", \editorname{}" * }
      if$
    }
  if$
}

FUNCTION {format.title}
{ title empty$
    { "" }
    % { title "t" change.case$ }
	{ title } % Je trouve que c'est mieux de laisser le titre tel quel, RK.
  if$
}

FUNCTION {n.dashify}
{ 't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
        { t #1 #2 substring$ "--" = not
%            { "--" *
            { "-" *
              t #2 global.max$ substring$ 't :=
            }
            {   { t #1 #1 substring$ "-" = }
                { "-" *
                  t #2 global.max$ substring$ 't :=
                }
              while$
            }
          if$
        }
        { t #1 #1 substring$ *
          t #2 global.max$ substring$ 't :=
        }
      if$
    }
  while$
}

FUNCTION {format.date}
{ year empty$
    { month empty$
        { "No year of publication in " cite$ * warning$ } %@@
        { "there's a month but no year in " cite$ * warning$
          month
        }
      if$
    }
    { month empty$
        'year
        { month " " * year * }
      if$
    }
  if$
}

FUNCTION {format.confdate} % date of conference
% @@@ A REVOIR
{ date empty$
    {"" }
    {date }
  if$
}

FUNCTION {format.btitle}
{ title emphasize
}

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$ * *
}

FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}

FUNCTION {format.bvolume}
{ volume empty$
    { "" }
    { "\volumename{}" volume tie.or.space.connect
      series empty$
        'skip$
        { " \ofname{} " * series emphasize * }
      if$
      "\volumename{} and \numbername{}" number either.or.check
    }
  if$
}

FUNCTION {format.number.series}
{ volume empty$
    { number empty$
        { series field.or.null }
        { output.state mid.sentence =
            { "\numbername{}" }
            { "\Numbername{}" }
          if$
          number tie.or.space.connect
          series empty$
            { "there's a number but no series in " cite$ * warning$ }
            { "\inname{}" * series * }
          if$
        }
      if$
    }
    { "" }
  if$
}

FUNCTION {format.edition}
{ edition empty$
    { "" }
    { output.state mid.sentence =
        { edition "l" change.case$ " \editionname{}" * }
        { edition "t" change.case$ " \editionname{}" * }
      if$
    }
  if$
}

INTEGERS { multiresult }

FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
        { #1 'multiresult := }
        { t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}

FUNCTION {format.pages}
{ pages empty$
    { "" }
    { pages multi.page.check
        { " \pagesname" pages n.dashify tie.or.space.connect }
        { " \pagename" pages tie.or.space.connect }
      if$
    }
  if$
}

FUNCTION {format.vol.num.year.pages} %%%%%%%%% NON UTILISÉ ICI
% utilisé exclusivemment par le type d'ouvrage = article
{ volume empty$
     { "there's no volume in " cite$ * warning$
       ""
     }
     { "\volumename{}" volume *}
  if$
  number empty$
    'skip$
    { ", " * "\numbername{}" number * *
      volume empty$
        { "there's a number but no volume in " cite$ * warning$ }
        'skip$
      if$
    }
  if$
  year empty$
     { " empty year in " cite$ * warning$ }
     {", " * year *}
  if$
  pages empty$
    'skip$
    { duplicate$ empty$
        { pop$ format.pages }
        { ", " * format.pages * }
      if$
    }
  if$
}

FUNCTION {format.vol.num.pages}
% utilisé exclusivemment par le type d'ouvrage = article
{ volume empty$
     { "there's no volume in " cite$ * warning$
       ""
     }
     { "\volumename{}" volume *}
  if$
  number empty$
    'skip$
    { ", " * "\numbername{}" number * *
      volume empty$
        { "there's a number but no volume in " cite$ * warning$ }
        'skip$
      if$
    }
  if$
  pages empty$
    'skip$
    { duplicate$ empty$
        { pop$ format.pages }
        { ", " * format.pages * }
      if$
    }
  if$
}

FUNCTION {format.chapter.pages}
{ chapter empty$
    'format.pages
    { type empty$
        { "\chaptername{}" }
        { type "l" change.case$ }
      if$
      chapter tie.or.space.connect
      pages empty$
        'skip$
        { ", " * format.pages * }
      if$
    }
  if$
}

FUNCTION {format.in.ed.booktitle}
{ booktitle empty$
    { "" }
    { editor empty$
        { "\Inname{}" booktitle emphasize * }
        { "\Inname{}" format.editors * ", " * booktitle emphasize * }
      if$
    }
  if$
}

FUNCTION {empty.misc.check}
{ author empty$ title empty$ howpublished empty$
  month empty$ year empty$ note empty$
  and and and and and
  key empty$ not and
    { "all relevant fields are empty in " cite$ * warning$ }
    'skip$
  if$
}

FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
%       type "t" change.case$
      type
    }
  if$
}

FUNCTION {format.tr.number}
{ type empty$
    { "\technicalreportname{}" }
    'type
  if$
  number empty$
    { "t" change.case$ }
    { "\numbername{}" number * tie.or.space.connect }
  if$
}

FUNCTION {format.article.crossref}
{ key empty$
    { journal empty$
        { "need key or journal for " cite$ * " to crossref " * crossref *
          warning$
          ""
        }
        { "\Inname{}\textit{" journal * "}" * }
      if$
    }
    { "\Inname{}" key * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {format.crossref.editor}
{ editor #1 "{vv~}{ll}" format.name$
  editor num.names$ duplicate$
  #2 >
    { pop$ " \etalter{}" * }
    { #2 <
        'skip$
        { editor #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
            { "\etalter{}" * }
            { "\andname{}" * editor #2 "{vv~}{ll}" format.name$ * }
          if$
        }
      if$
    }
  if$
}

FUNCTION {format.book.crossref}
{ volume empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      "\Inname{}"
    }
    { "\Volumename{}" volume tie.or.space.connect
      " of " *
    }
  if$
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
        { series empty$
            { "need editor, key, or series for " cite$ * " to crossref " *
              crossref * warning$
              "" *
            }
            { "\textit{" * series * "}" * }
          if$
        }
        { key * }
      if$
    }
    { format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {format.incoll.inproc.crossref}
{ editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
        { booktitle empty$
            { "need editor, key, or booktitle for " cite$ * " to crossref " *
              crossref * warning$
              ""
            }
            { "\Inname{}\textit{" booktitle * "}" * }
          if$
        }
        { "\Inname{}" key * }
      if$
    }
    { "\Inname{}" format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}

% Types de references %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTION {ARTICLE}
{ output.bibitem
  format.authors "author" output.check
  format.title "title" output.check.guillemet
  new.block
  crossref missing$
    { journal emphasize "journal" output.check
%      format.vol.num.year.pages output
      format.vol.num.pages output
    }
    { format.article.crossref output.nonnull
      format.pages output
    }
  if$
  publisher output
  format.date "year" output.check
  new.block
  note output
  fin.entry
}

FUNCTION {BOOK}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
    }
  if$
  format.btitle "title" output.check
  crossref missing$
    { format.bvolume output
      new.block
      format.number.series output
      new.sentence
      publisher "publisher" output.check
      address output
    }
    { new.block
      format.book.crossref output.nonnull
    }
  if$
  format.edition output
  format.date "year" output.check
  new.block
  note output
  fin.entry
}

FUNCTION {BOOKLET}
{ output.bibitem
  format.authors output
  format.title "title" output.check
  howpublished address new.block.checkb
  howpublished output
  address output
  format.date output
  new.block
  note output
  fin.entry
}

FUNCTION {INBOOK}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
    }
  if$
  format.btitle "title" output.check.guillemet
  crossref missing$
    { format.bvolume output
      format.chapter.pages "chapter and pages" output.check
      new.block
      format.number.series output
      new.sentence
      publisher "publisher" output.check
      address output
    }
    { format.chapter.pages "chapter and pages" output.check
      new.block
      format.book.crossref output.nonnull
    }
  if$
  format.edition output
  format.date "year" output.check
  new.block
  note output
  fin.entry
}

FUNCTION {INCOLLECTION}
{ output.bibitem
  format.authors "author" output.check
  format.title "title" output.check.guillemet
  new.block
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.bvolume output
      format.number.series output
      format.chapter.pages output
      new.sentence
      publisher "publisher" output.check
      address output
      format.edition output
      format.date "year" output.check
    }
    { format.incoll.inproc.crossref output.nonnull
      format.chapter.pages output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {INPROCEEDINGS}
{ output.bibitem
  format.authors "author" output.check
  format.title "title" output.check.guillemet
  new.block
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.bvolume output
      format.number.series output
      address empty$
        { organization publisher new.sentence.checkb
          organization output
          publisher output
%          format.confdate "date of conference" output.check 
          date output
        }
        { address output.nonnull
%          format.confdate "date of conference" output.check
          date output
          new.sentence
          organization output
          publisher output
        }
      if$
      format.pages output
    }
    { format.incoll.inproc.crossref output.nonnull
      format.pages output
    }
  if$
  format.date "year" output.check
  new.block
  note output
  fin.entry
}

FUNCTION {CONFERENCE} {INPROCEEDINGS}

FUNCTION {MANUAL}
{ output.bibitem
  author empty$
    { organization empty$
        'skip$
        { organization output.nonnull
          address output
        }
      if$
    }
    { format.authors output.nonnull }
  if$
  format.title "title" output.check %.guillemet
  author empty$
    { organization empty$
        { address new.block.checka
          address output
        }
        'skip$
      if$
    }
    { organization address new.block.checkb
      organization output
      address output
    }
  if$
  format.edition output
  format.date output
  new.block
  note output
  fin.entry
}

FUNCTION {MASTERSTHESIS}
{ output.bibitem
  format.authors "author" output.check
  format.title "title" output.check %.guillemet
  new.block
  "Master's thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  format.date "year" output.check
  new.block
  note output
  fin.entry
}

FUNCTION {MISC}
{ output.bibitem
  format.authors output
  title howpublished new.block.checkb
  format.title output.guillemet
  howpublished new.block.checka
  howpublished output
  format.date output
  new.block
  note output
  fin.entry
  empty.misc.check
}

FUNCTION {PHDTHESIS}
{ output.bibitem
  format.authors "author" output.check
  format.title "title" output.check% .guillemet
  new.block
  "PhD thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  format.date "year" output.check
  new.block
  note output
  fin.entry
}

FUNCTION {PROCEEDINGS}
{ output.bibitem
  editor empty$
    { organization output }
    { format.editors output.nonnull }
  if$
  format.btitle "title" output.check
  format.bvolume output
  format.number.series output
  address empty$
    { editor empty$
        { publisher new.sentence.checka }
        { organization publisher new.sentence.checkb
          organization output
        }
      if$
      publisher output
%      format.confdate "date of conference" output.check
      date output
    }
    { address output.nonnull
%      format.confdate "date of conference" output.check
      date output
      new.sentence
      editor empty$
        'skip$
        { organization output }
      if$
      publisher output
    }
  if$
  format.date "year" output.check
  new.block
  note output
  fin.entry
}

FUNCTION {TECHREPORT}
{ output.bibitem
  format.authors "author" output.check
  format.title "title" output.check %.guillemet
  new.block
  format.tr.number output.nonnull
  institution "institution" output.check
  address output
  format.date "year" output.check
  new.block
  note output
  fin.entry
}

FUNCTION {UNPUBLISHED}
{ output.bibitem
  format.authors "author" output.check
  format.title "title" output.check %.guillemet
  new.block
  note "note" output.check
  format.date output
  fin.entry
}

FUNCTION {default.type} { misc }

% MACROS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MACRO {jan} {"\Jan{}"}
MACRO {feb} {"\Feb{}"}
MACRO {mar} {"\Mar{}"}
MACRO {apr} {"\Apr{}"}
MACRO {may} {"\May{}"}
MACRO {jun} {"\Jun{}"}
MACRO {jul} {"\Jul{}"}
MACRO {aug} {"\Aug{}"}
MACRO {sep} {"\Sep{}"}
MACRO {oct} {"\Oct{}"}
MACRO {nov} {"\Nov{}"}
MACRO {dec} {"\Dec{}"}

% Actions standards %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

READ

FUNCTION {sortify}
{ purify$
  "l" change.case$
}

INTEGERS { len }

FUNCTION {chop.word}
{ 's :=
  'len :=
  s #1 len substring$ =
    { s len #1 + global.max$ substring$ }
    's
  if$
}

INTEGERS { et.al.char.used }

FUNCTION {initialize.et.al.char.used}
{ #0 'et.al.char.used :=
}

EXECUTE {initialize.et.al.char.used}

% On prend les 3 premieres lettres du premier auteur en majuscule comme cle:

FUNCTION {format.lab.names}
{ 's :=
	s #1 "{ll{}}" format.name$	% duplicate$ top$
	#3 text.prefix$			% duplicate$ top$
	"u" change.case$		% duplicate$ top$
	"~" *				% duplicate$ top$
}

FUNCTION {author.key.label}
{ author empty$
    { key empty$
        { cite$ #1 #3 substring$ }
        { key #3 text.prefix$ }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {author.editor.key.label}
{ author empty$
    { editor empty$
        { key empty$
            { cite$ #1 #3 substring$ }
            { key #3 text.prefix$ }
          if$
        }
        { editor format.lab.names }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {author.key.organization.label}
{ author empty$
    { key empty$
        { organization empty$
            { cite$ #1 #3 substring$ }
            { "The " #4 organization chop.word #3 text.prefix$ }
          if$
        }
        { key #3 text.prefix$ }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {editor.key.organization.label}
{ editor empty$
    { key empty$
        { organization empty$
            { cite$ #1 #3 substring$ }
            { "The " #4 organization chop.word #3 text.prefix$ }
          if$
        }
        { key #3 text.prefix$ }
      if$
    }
    { editor format.lab.names }
  if$
}

FUNCTION {calc.label}
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.key.label
    { type$ "proceedings" =
        'editor.key.organization.label
        { type$ "manual" =
            'author.key.organization.label
            'author.key.label
          if$
        }
      if$
    }
  if$
  duplicate$
  year field.or.null purify$ #-1 #2 substring$
  *
  'label :=
  year field.or.null purify$ #-1 #4 substring$
  *
  sortify 'sort.label :=
}

FUNCTION {sort.format.names}
{ 's :=
  #1 'nameptr :=
  ""
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { nameptr #1 >
        { "   " * }
        'skip$
      if$
      s nameptr "{vv{ } }{ll{ }}{  ff{ }}{  jj{ }}" format.name$ 't :=
      nameptr numnames = t "others" = and
        { "\etalter{}" * }
        { t sortify * }
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {sort.format.title}
{ 't :=
  "A " #2
    "An " #3
      "The " #4 t chop.word
    chop.word
  chop.word
  sortify
  #1 global.max$ substring$
}

FUNCTION {author.sort}
{ author empty$
    { key empty$
        { "to sort, need author or key in " cite$ * warning$
          ""
        }
        { key sortify }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {author.editor.sort}
{ author empty$
    { editor empty$
        { key empty$
            { "to sort, need author, editor, or key in " cite$ * warning$
              ""
            }
            { key sortify }
          if$
        }
        { editor sort.format.names }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {author.organization.sort}
{ author empty$
    { organization empty$
        { key empty$
            { "to sort, need author, organization, or key in " cite$ * warning$
              ""
            }
            { key sortify }
          if$
        }
        { "The " #4 organization chop.word sortify }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {editor.organization.sort}
{ editor empty$
    { organization empty$
        { key empty$
            { "to sort, need editor, organization, or key in " cite$ * warning$
              ""
            }
            { key sortify }
          if$
        }
        { "The " #4 organization chop.word sortify }
      if$
    }
    { editor sort.format.names }
  if$
}

FUNCTION {presort}
{ calc.label
  sort.label
  "    "
  *
  type$ "book" =
  type$ "inbook" =
  or
    'author.editor.sort
    { type$ "proceedings" =
        'editor.organization.sort
        { type$ "manual" =
            'author.organization.sort
            'author.sort
          if$
        }
      if$
    }
  if$
  *
  "    "
  *
  year field.or.null sortify
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {presort}

SORT

STRINGS { longest.label last.sort.label next.extra }

INTEGERS { longest.label.width last.extra.num }

FUNCTION {initialize.longest.label}
{ "" 'longest.label :=
  #0 int.to.chr$ 'last.sort.label :=
  "" 'next.extra :=
  #0 'longest.label.width :=
  #0 'last.extra.num :=
}

FUNCTION {forward.pass}
{ last.sort.label sort.label =
    { last.extra.num #1 + 'last.extra.num :=
      last.extra.num int.to.chr$ 'extra.label :=
    }
    { "a" chr.to.int$ 'last.extra.num :=
      "" 'extra.label :=
      sort.label 'last.sort.label :=
    }
  if$
}

FUNCTION {reverse.pass}
{ next.extra "b" =
    { "a" 'extra.label := }
    'skip$
  if$
  label extra.label * 'label :=
  label width$ longest.label.width >
    { label 'longest.label :=
      label width$ 'longest.label.width :=
    }
    'skip$
  if$
  extra.label 'next.extra :=
}

EXECUTE {initialize.longest.label}

ITERATE {forward.pass}

REVERSE {reverse.pass}

FUNCTION {begin.bib}
{ et.al.char.used
    { "\newcommand{\etalchar}[1]{$^{#1}$}" write$ newline$ }
    'skip$
  if$
  preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{}" write$ newline$
}

EXECUTE {begin.bib}

EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}

EXECUTE {end.bib}
