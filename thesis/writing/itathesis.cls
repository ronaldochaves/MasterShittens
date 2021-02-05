% !TeX encoding = ISO-8859-1

%==============================================================================
% This is the ITA LateX class - classe de documento com padrões do ITA
% written by
%              Fábio Fagundes Silveira   -  ffs [at] ita [dot] br
%
%              Benedito C. O. Maciel     -  bcmaciel [at] ita [dot] br
%              Giovani Volnei Meinertz   -  giovani [at] ita [dot] br
%              Tarcisio A. B. Gripp      -  tarcisio.gripp [at] gmail [dot] com
%
% last updated: 
%        $Id: ita.cls 42 2015-04-08$
%        $HeadURL: http://www.apgita.org.br/apgita/teses-e-latex.php $
%
%
% Copyright (C) 2004-2008 ITALUS - ITA LaTeX Users Group "Society"  :-)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% ITALUS
% Instituto Tecnologico de Aeronautica --- ITA
% Sao Jose dos Campos, Brazil
% HomePages:        http://www.comp.ita.br/italus
%                   http://groups.yahoo.com/group/italus/
% Discussion list: italus {at} yahoogroups.com
%==============================================================================

%==============================================================================
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{itathesis}[2006/03/10 v2.0 Classe de documento com padrões do ITA]

% % --------------------------------------------------------------------
% %  Verificação se alvo será gerado pelo PDFLaTeX ou não
% % --------------------------------------------------------------------
\IfFileExists{ifpdf.sty}{
	\RequirePackage{ifpdf}
}{\newif\ifpdf}
%\newif\ifpdf
%  \ifx\pdfoutput\undefined
 %     \pdffalse
  %\else
   %   \pdfoutput=1
    %  \pdftrue
 %\fi

%==============================================================================
%%% Declaração das opções
%==============================================================================

\DeclareOption{defaults}{
	\input{defaults.def}
}

% Graduação (TG)
\DeclareOption{tg}{
	\input{tg.def}
}

% Mestrado
\DeclareOption{msc}{
	\input{msc.def}
}

% Doutorado
\DeclareOption{dsc}{
	\input{dsc.def}
}

% Exame de Qualificação
\DeclareOption{quali}{
	\input{quali.def}
}

% Versão em Inglês
\newif\if@eng
\DeclareOption{eng}{
    \@engtrue
}

\DeclareOption{classic}{
	\newcommand{\itaformatworktitle}[1]{\Large \textbf{#1}}
	\newcommand{\itaformatauthorname}[1]{\large \textbf{#1}}
	\newcommand{\itaformattitletext}[1]{\large{#1}}
	\newcommand{\itaformatadvisorword}[1]{#1}
	\renewcommand{\itaspacing}{1.5}
}

%==============================================================================
%%% draft: 'draft version': imprime ´Versão Preliminar´ + data da impressão
%==============================================================================
\newif\if@draft
\DeclareOption{draft}{
    \@drafttrue
}

\ExecuteOptions{defaults, classic}

\ProcessOptions

\AtEndOfClass{
    \renewcommand{\baselinestretch}{\itaspacing}
}
  
%==============================================================================
%%% Herança e Pacotes da classe book.cls
%==============================================================================
\LoadClass[12pt, a4paper, oneside]{book}
\RequirePackage{calc}
\RequirePackage{ifthen}

\if@eng
	\RequirePackage[english]{babel}
	\input{itaeng.babel}
\else
	% Comandos em português (e.g. faz \chapter virar Capítulo no texto) e hifenização correta
	\RequirePackage[english, portuguese]{babel}
	\input{itabr.babel}
\fi

\RequirePackage[latin1]{inputenc}		% Muda a condificação da fonte, possibilitando copiar textos com acentos etc
\RequirePackage{setspace}
\RequirePackage{fancyhdr}
\RequirePackage{amsfonts}               % Fontes p/ conjuntos: R, Z, H etc
\RequirePackage{lastpage}
\RequirePackage{indentfirst}
%ajustando fonte das legendas. Se quiser aumentar um pouco as legendas, troque footnotesize por small
%\RequirePackage[labelsep = endash, font = footnotesize]{caption}
\RequirePackage[labelsep = endash]{caption}
\RequirePackage{longtable}
\RequirePackage{array}
\RequirePackage{tabularx}
\RequirePackage[en-US, calc]{datetime2}

%==============================================================================
%%% Glossário
\usepackage{makeglo}
\renewcommand{\glossmaxlabelwidth}{\hspace{3cm}}

%% Opções de cores
%% linkcolor = blue, citecolor = gray, urlcolor = cyan
\RequirePackage{color}        % Aplicar cores no texto
\RequirePackage{hyphenat}     % Hifenização
\RequirePackage{makeidx}      % Índice

%% Hyperref
%==============================================================================
\definecolor{gray}{rgb}{0.4, 0.4, 0.4}
\definecolor{bblue}{rgb}{0, 0, 0.8} 		% Tom azul mais escuro

%\RequirePackage[debug, a4paper, pagebackref = true, breaklinks, bookmarksopen = true, colorlinks = true, linkcolor = bblue, citecolor = bblue, urlcolor = bblue, bookmarksnumbered = true, hyperindex = true, verbose = true]{hyperref}
% Gera pdf/dvi com links

% Links alterados para cor preta
\RequirePackage[debug, a4paper, pagebackref = true, breaklinks, bookmarksopen = true, colorlinks = true, linkcolor = black, citecolor = black, urlcolor = black, bookmarksnumbered = true, hyperindex = true, verbose = true]{hyperref} 
%==============================================================================

\ifpdf
	\RequirePackage{ae} 					% Fontes para pdf
	\RequirePackage[pdftex]{graphicx}
	\RequirePackage[pdftex]{geometry}
\else
	\RequirePackage[T1]{fontenc}			% Possibilita usar "á" ao invés de "\'a", "ç" ao invés de "\c{c}" etc
	\RequirePackage[dvips]{graphicx}
	\RequirePackage[dvips]{geometry}
\fi

% Dimensões das páginas
\geometry{a4paper, left = 3cm, right = 2cm, top = 3cm, bottom = 2cm, verbose, headheight = 0.56cm, headsep = 0.5cm, footskip = 3em}

\newif\ifhyperref
\@ifpackageloaded{hyperref}{\hyperreftrue}{\hyperreffalse}

%==============================================================================
%%% ABNTeX
\RequirePackage[alf, abnt-emphasize = bf, abnt-etal-cite = 2, abnt-etal-list = 0, abnt-etal-text = it]{abntcite}
%\RequirePackage{natbib}

%==============================================================================
% ALTERAR TAMANHO DA FONTE DOS TÍTULOS
%\RequirePackage{titlesec} 							% Alterar tamanho da fonte
%\titleformat{\section}
% {\large\bfseries}{\thesection}{0.5em}{}
%\titleformat{\subsection}
% {\normalsize\bfseries}{\thesubsection}{0.5em}{}	
%\titleformat{\subsubsection}
% {\normalsize\bfseries}{\thesubsubsection}{0.5em}{}	
%\titleformat{\paragraph} 							% para 4o nível
% {\normalsize\bfseries}{\theparagraph}{0.5em}{}	

%==============================================================================
% Redefinição da data (para ter dia, mês e ano separados)
%==============================================================================
\DTMnewdatestyle{monthdeyyyy}% label    % TODO add this to Portuguese Module for datetime2 Package
{% definitions
	\selectlanguage{portuguese}
	\renewcommand*{\DTMdisplaydate}[4]{%
		\DTMmonthname{##2}\DTMportugesmonthyearsep##1}%
	\renewcommand*{\DTMDisplaydate}{\DTMdisplaydate}%
}

\renewcommand{\date}[3]{
	\DTMsavedate{defDate}{#3-#2-#1}
	\gdef\@publday{#1}
	\gdef\@publmonth{#2}
	\gdef\@publyear{#3}
}

\newcommand{\FRDitadata}[3]{
	\DTMsavedate{publDate}{#3-#2-#1}
	\gdef\@FRDday{#1}
	\gdef\@FRDmonth{#2}
	\gdef\@FRDyear{#3}
}

%==============================================================================
%%% Novos comandos: \inst, \author, \authorAddress, \advisor, \coadvisor, \boss etc
%==============================================================================

% Instituição
\newcommand{\inst}[2]{               	 % Nome e abreviação
	\gdef\@instFull{#1}
	\gdef\@instAbbr{#2}
}

% Endereço da instituição
\newcommand{\instAddress}[5]{               	 % Rua, CEP, Cidade e Estado
	\gdef\@instStreet{#1}
	\gdef\@instZIP{#2}
	\gdef\@instCity{#3}
	\gdef\@instState{#4}
	\gdef\@instCampus{#5}
}

% Autor
\renewcommand{\author}[3]{               % Nome, Sobrenome e Gênero
	\gdef\@authorFirst{#1}
	\gdef\@authorLast{#2}
	\gdef\@authorGender{#3}
	\gdef\@authorRole{\@authorword}
}

% Endereço do autor (utilizado no verso da folha de rosto)
\renewcommand{\authorAddress}[4]{        % Rua, CEP, Cidade e Estado
	\gdef\@authorStreet{#1}
	\gdef\@authorZIP{#2}
	\gdef\@authorCity{#3}
	\gdef\@authorState{#4}
}

% Orientador
\newcommand{\advisor}[5]{				 % Profissão, Título, Nome, Gênero e Instituição
	\gdef\advisorCarrer{#1}
	\gdef\advisorDegree{#2}
	\renewcommand{\advisorName}{#3}
	\gdef\advisorGender{#4}
	\gdef\advisorInst{#5}
	\gdef\advisorRole{\itaadvisorword}
	\renewcommand{\advisorTitle}{\genderCarrer{\advisorCarrer}{\advisorGender}~\genderDegree{\advisorDegree}{\advisorGender}}
	\gdef\@advisor{\advisorTitle~\advisorName \\ \advisorRole}
	\gdef\@advisorTG{\advisorRole \\ \advisorTitle~\advisorName\ (\advisorInst)}
}

% Coorientador (opcional)
\newcommand{\coadvisor}[5]{				 % Profissão, Título, Nome, Gênero e Instituição
	\gdef\coadvisorCarrer{#1}
	\gdef\coadvisorDegree{#2}
	\renewcommand{\coadvisorName}{#3}
	\gdef\coadvisorGender{#4}
	\gdef\coadvisorInst{#5}
	\gdef\coadvisorRole{\itacoadvisorword}
	\renewcommand{\coadvisorTitle}{\genderCarrer{\coadvisorCarrer}{\coadvisorGender}~\genderDegree{\coadvisorDegree}{\coadvisorGender}}
	\gdef\@coadvisor{\coadvisorTitle~\coadvisorName \\ \coadvisorRole}
	\gdef\@coadvisorTG{\coadvisorRole \\ \coadvisorTitle~\coadvisorName\ (\coadvisorInst)}
}

% Pró-Reitor de Pós-Graduação e Pesquisa
\newcommand{\boss}[4]{					 % Profissão, Título, Nome e Gênero
	\gdef\bossCarrer{#1}
	\gdef\bossDegree{#2}
	\renewcommand{\bossName}{#3}
	\gdef\bossGender{#4}
	\gdef\bossRole{\ifthenelse{\equal{\itaoptionlevel}{tg}}{\itabosscourseword}{\itabossword}}
	\renewcommand{\bossTitle}{\genderCarrer{\bossCarrer}{\bossGender}~\genderDegree{\bossDegree}{\bossGender}}
	\gdef\@boss{\bossTitle~\bossName \\ \bossRole}
}

% Curso de Graduação/Programa de Pós-Graduação
\newcommand{\course}[1]{
	\renewcommand{\itaworkcourse}{#1}
}

% Área do Curso													    

% Coorientador
\newcommand{\area}[1]{
	\renewcommand{\itaworkarea}{#1}
}

% Número de CDU
\newcommand{\cdu}[1]{
	\renewcommand{\cdunumber}{#1}
}

%==============================================================================
%%% Hack: Página de Rosto para Teses
%==============================================================================
\newcommand{\itatitletext}{\itatitletextbabel}

\renewcommand{\maketitle}{%
   \begin{titlepage}
      \ifthenelse{\equal{\itaoptionlevel}{tg}}{\maketgcover}{}
      \pdfbookmark[0]{\nametitlepage}{titulo}
      
      % verifica se é uma tese ou um exame de qualificação
      \ifthenelse{\equal{\itaoptionlevel}{quali}}{\thequalititlepage}{\ifthenelse{\equal{\itaoptionlevel}{tg}}{\thetgtitlepage}{\thethesistitlepage}}
   \end{titlepage}

   \ifthenelse{\equal{\itaoptionlevel}{quali}}{\setcounter{page}{2}}{\setcounter{page}{4}}	
}

%==============================================================================
% Folha de Rosto do Exame de Qualificação para Doutorado
%==============================================================================
\newcommand{\thequalititlepage}{
	
	\noindent
	\center
	\begin{spacing}{1.5}
		{\Large \itanametitlequalibabel } 		\\		%  \@instAbbr -- \textsc{\@instFull}} \\
		{\Large \itaposnametitlequalibabel }	\\		%  \textsc{Pós-Graduação em \itaworkcourse}} \\
		{\large \itaarenametitlequalibabel }			%  \textsc{Área de \itaworkarea}}
	\end{spacing}

	\vfill

	\begin{spacing}{1.2}
		\itaformatworktitle{\MakeUppercase{\@title}}
	\end{spacing}

	\vfill

	\if@eng
		by\\
	\else
		por\\
	\fi
	\itaformatauthorname{\@authorFirst\ \@authorLast}

	\vfill

	\normalsize
	\hspace{\fill}
	\parbox[b]{8.5cm}{\itaqualifictext}

	\vfill

	\begin{spacing}{1.5}
		\@advisor
	\end{spacing}

	\@ifundefined{@coadvisor}{}{      					% Co-advisor é opcional
		\vfill
		\begin{spacing}{1.5}
			\@coadvisor
		\end{spacing}
	}

	\vfill

	\if@eng
		\begin{spacing}{1.5}
			\@instCity: \DTMusedate{defDate}
		\end{spacing}
	\else
		\begin{spacing}{1.5}
			\DTMsetdatestyle{monthdeyyyy}
			\@instCity, \DTMusedate{defDate}
		\end{spacing}
	\fi
}

%==============================================================================
% Folha de Rosto do Trabalho de Graduação
%==============================================================================
\newcommand{\thetgtitlepage}{
	
	\newpage
	\thispagestyle{empty}  						% Suprime o cabeçalho e numeração da página
	\noindent
	
	\begin{flushright}							% Número de CDU
		\MakeUppercase{\footnotesize CDU \cdunumber}
	\end{flushright}

	\center

	\itaformatauthorname{\@authorFirst\ \@authorLast}
	\vfill

	\vspace*{\stretch{1}}

	\begin{spacing}{1.0}
		\itaformatworktitle{\MakeUppercase{\@title}}
	\end{spacing}

	\vspace*{\stretch{2}}

	\begin{spacing}{1.5}
		\@advisorTG
	\end{spacing}

	\@ifundefined{@coadvisor}{}{      			% Co-advisor é opcional
		\vspace*{\stretch{2}}
		\begin{spacing}{1.5}
			\@coadvisorTG
		\end{spacing}
	}

	\vspace*{\stretch{3}}

	\itaformatauthorname{\MakeUppercase{\textbf{\itaworkcourse}}}

	\vspace*{\stretch{2}}

	\vspace*{\stretch{2}}

	\begin{spacing}{1.0}
		\textsc{\Large{S}\large{ão} \Large{J}\large{osé} dos \Large{C}\large{ampos}} \\ \textsc{\Large{i}\large{nstituto}{} \Large{t}\large{ecnológico}{} de \Large{a}\large{eronáutica}{}}
		
		\vspace*{\stretch{2}}
		
		\DTMfetchyear{defDate}
	\end{spacing}

	\@makecippage        						% Cria verso da folha de rosto
	\@listexaminators    						% Cria página de membros da banca
}

%==============================================================================
% Folha de Rosto das Teses de Mestrado e Doutorado
%==============================================================================
\newcommand{\thethesistitlepage}{
	
	\noindent

	\begin{spacing}{1.2}
		\noindent \parbox[c]{\textwidth}{\itaformattitletext{\itatitletext}}
	\end{spacing}

	\vspace*{\stretch{2}}

	\center

	\itaformatauthorname{\@authorFirst\ \@authorLast}
	
	\vfill
	\vspace*{\stretch{1}}

	\begin{spacing}{1.2}
		\itaformatworktitle{\MakeUppercase{\@title}}
	\end{spacing}

	\vspace*{\stretch{2}}

	\normalsize
	\itathesisapprovedmembersbabel

	\vspace*{\stretch{3}}

	\begin{spacing}{1.5}
		\@advisor
	\end{spacing}

	\@ifundefined{@coadvisor}{}{      			% Co-advisor é opcional
		\vspace*{\stretch{2}}
		\begin{spacing}{1.5}
			\@coadvisor
		\end{spacing}
	}

	\vspace*{\stretch{2}}

	\begin{spacing}{1.5}
		\@boss
	\end{spacing}

	\vspace*{\stretch{2}}

	\begin{spacing}{1.0}
		\@instCampus \\ \@instCity, \@instState\ -- \@instCountry \\ \DTMfetchyear{defDate}
	\end{spacing}

	\@makecippage        						% Cria verso da folha de rosto
	\@listexaminators    						% Cria página de membros da banca
}

%==============================================================================
% Capa do TG
%==============================================================================

\newcommand{\maketgcover}{
	
	\pdfbookmark[0]{\namecover}{cover}
	\thispagestyle{empty}

	\hspace{-0.8cm}\fbox{
		\parbox[c][660pt][c]{\textwidth}{
			\centering        
		
			\itaformatworktitle{\MakeUppercase{\@instFull}}
			
			\vspace*{\stretch{2}}
			
			\includegraphics[width = 0.5\textwidth]{Cap0/ita_logo}
			
			\vspace*{\stretch{2}}
			
			\itaformatauthorname{\@authorFirst\ \@authorLast}
			
			\vspace*{\stretch{1}}
			
			\begin{spacing}{1.0}
				\itaformatworktitle{\MakeUppercase{\@title}}
			\end{spacing}
			
			\vspace*{\stretch{3}}
			
			\begin{spacing}{1.0}
				\Large{\finalpaper \\
				\DTMfetchyear{defDate}}
			\end{spacing}
			
			\vspace*{\stretch{2}}
			
			\Huge{\textbf{\courseof\ \itaworkcourse}}
			
			\vspace*{\stretch{1}}
			\vfill
		}
	}
}

%==============================================================================
% CIP - Verso da Folha de Rosto
%==============================================================================
\newcommand{\@makecippage}{
	
	\newpage
	\pdfbookmark[0]{\namecip}{cip}
	\thispagestyle{empty}
	
	\begin{flushleft}
		\begin{footnotesize}
			\begin{onehalfspacing}
				{
					\setlength{\baselineskip}{0.8\baselineskip}
					\textbf{\itaciptitlebabel} 		\\
					\textbf{\itacipbibbabel}		\\
				}		
				
				\@cipwindow
				\@bibref
				\@rights
			\end{onehalfspacing}
		\end{footnotesize}
	\end{flushleft}
	
	\vfill
}

% Quadro com informações CIP
\def\@cipwindow{
	\framebox[\textwidth]{
		\begin{minipage}{\textwidth - 10mm}
			\begin{flushleft}
				\vspace*{0.4ex}
				\setlength{\parindent}{0.5em}
				
				\begin{footnotesize}
					\begin{onehalfspacing}
						\setlength{\baselineskip}{0.7\baselineskip}
						\noindent\@authorLast, \@authorFirst\par
						\@title\ / \@authorFirst\ \@authorLast.			\\
						\noindent \@instCity, \DTMfetchyear{defDate}.\par
						\NumPenultimaPagina f.\par
						\vspace{\baselineskip}
						\ifthenelse{\equal{\itaoptionlevel}{tg}}{\itaciptgdata1babel}{\itacipdata1babel}
						\vspace{0.4ex}
					\end{onehalfspacing}
				\end{footnotesize}
			\end{flushleft}
		\end{minipage}
	}
}

% Referência Bibliográfica
\def\@bibref{%
	
	\vspace*{\stretch{2}}
	\normalsize
	\textbf{\itarefbibnamecipbabel} \par
	
	\begin{onehalfspacing}
		{
			\setlength{\baselineskip}{0.7\baselineskip}
			\MakeUppercase{\@authorLast}, \@authorFirst. \textbf{\@title}. \DTMfetchyear{defDate}.
			\NumPenultimaPagina f. \ifthenelse{\equal{\itaoptionlevel}{tg}}{\itaciptg2babel}{\itacip2babel} \\
		}
	\end{onehalfspacing}
}

% Cessão de Direitos
\def\@rights{%
	\begin{onehalfspacing}
		
		\vspace*{\stretch{4}}
		
		{\normalsize\textbf{\itacedirnamebabel}} \par
		{
			\setlength{\baselineskip}{0.7\baselineskip}
			{\footnotesize \itacipauthorname} {\normalsize \@authorFirst\ \@authorLast} \\
			{\footnotesize \itacipworktitlebabel} {\normalsize \@title.} 				\\
			{\footnotesize \itacipworktype} {\normalsize\ifthenelse{\equal{\itaoptionlevel}{tg}}{\itatglevel}{\itaworklevel} / \DTMfetchyear{defDate}}
			\\
		}
		
		\vspace{\baselineskip}
		
		\parbox[c]{\textwidth}{
			\normalsize\setlength{\baselineskip}{0.7\baselineskip} \ifthenelse{\equal{\itaoptionlevel}{tg}}{\itacesdirtgtextbabel}{\itacesdirtextbabel} \\
		}
		
		\vspace*{\stretch{2}}
		\rule{75mm}{0.2mm}\\
		{\normalsize\setlength{\baselineskip}{0.7\baselineskip}
		\@authorFirst\ \@authorLast\\
		\@authorStreet\\
		\@authorZIP\ -- \@authorCity\ -- \@authorState\\}
	\end{onehalfspacing}
}

%==============================================================================
%%% Banca Examinadora
%==============================================================================
\newcommand{\@listexaminators}
{
   \newpage
   \pdfbookmark[0]{\namepageexaminersbookmark}{aprovacao}
   \thispagestyle{empty}

   \noindent
   \ifthenelse{\equal{\itaoptionlevel}{tg}}
   {
   		\begin{spacing}{1.0}
   			\begin{center}
				\itaformatworktitle{\MakeUppercase{\@title}}
				\vfill
				\small{\itathesisapprovedmemberstg}
				\vspace*{\stretch{3}}
				
				\hspace{3cm}\hrulefill\hspace{3cm}
				\par 
				\small{\@authorFirst\ \@authorLast
				\par
				\@authorRole}
				\vspace*{\stretch{1}}
				 
				\hspace{3cm}\hrulefill\hspace{3cm}
				\par 
				\small{\advisorTitle~\advisorName\ (\advisorInst)
				\par
				\advisorRole}
				\vspace*{\stretch{1}}
				 		
				\@ifundefined{@coadvisor}{}{      % Co-advisor é opcional
				\hspace{3cm}\hrulefill\hspace{3cm}
				\par 
				\small{\coadvisorTitle~\coadvisorName\ (\coadvisorInst)
				\par
				\coadvisorRole}
				\vspace*{\stretch{1}}
				}
					
				\hspace{3cm}\hrulefill\hspace{3cm}
				\par 
				\small{\@boss}
				\vspace*{\stretch{2}}
					
				\if@eng
							\small{\@instCity: \DTMusedate{defDate}.}
				\else
							\small{\@instCity, {\selectlanguage{portuguese} \DTMusedate{defDate}}}
				\fi
			\end{center}
		\end{spacing}
			 
	} 
	{
		\begin{spacing}{1.0}
			\begin{center}
				\itaformatworktitle{\MakeUppercase{\@title}}
				\vfill
				\itaformatauthorname{\@authorFirst\ \@authorLast}
				\vfill
			\end{center}
			
			\raggedright{\normalsize{\namepageexaminers}}
			\vfill
			\raggedright\@memberslist
		\end{spacing}
		\vfill
		\center
		\large\textbf{\@instAbbr}
	}
}

%==============================================================================
%%% Dedicatória
%==============================================================================
\newenvironment{itadedication}{
	\thispagestyle{empty}%
	\pdfbookmark[0]{\namededicatory}{dedication}
	\begin{flushright}%
		\null\vfil \@beginparpenalty\@lowpenalty
		\vspace{120mm}
		\hspace{80mm}
		\begin{minipage}[b]{0.45\textwidth}
			\begin{spacing}{1.5}
}
{
			\end{spacing}
		\end{minipage}
	\end{flushright}\par\vfil\null\newpage
}

%==============================================================================
%%% Agradecimentos
%==============================================================================
%%alterado em 09/01/2015: retirado o alinhamento do texto à direta 
%\newenvironment{itathanks}[1][\nameacknowledgment]
%{
% \chapter*{#1}
% \pdfbookmark[0]{\nameacknowledgment}{acknowledgment}
% \begin{flushright}\begin{spacing}{1.5}}
% {\end{spacing}\end{flushright} \par\vfil\null\newpage}

\newenvironment{itathanks}[1][\nameacknowledgment]{
	\chapter*{#1}
	\pdfbookmark[0]{\nameacknowledgment}{acknowledgment}
}
{
	\par\vfil\null\newpage
}

%==============================================================================
%%% Epígrafe
%==============================================================================
\newenvironment{itaepigraph}{
	\thispagestyle{empty}%
	\ifhyperref\pdfbookmark[0]{\nameepigraph}{epigraph}\fi
	\begin{flushright}%
		\begin{spacing}{1}
			\mbox{}\vfill
			{\sffamily\itshape
}
{
			}
			\\ --- \textsc{Sir Isaac Newton}
		\end{spacing}
	\end{flushright}\par\vfil\null\newpage
}

%==============================================================================
%%% Resumo
%==============================================================================
\newenvironment{abstract}{
	\chapter*{\nameabstract}
 	\pdfbookmark[0]{\nameabstract}{resumo}
}
{
}

%==============================================================================
%%% Abstract
%==============================================================================
\newenvironment{englishabstract}{
	\chapter*{\nameengabstract}
	\pdfbookmark[0]{\nameengabstract}{abstract}
}
{
}

%==============================================================================
%%% Redefinição de chapter
%==============================================================================
\renewcommand\chapter{

	\if@openright
		\cleardoublepage
	\else
		\clearpage
	\fi
    
    \thispagestyle{empty}%
    \global\@topnum\z@
    \@afterindenttrue
    \secdef\@chapter\@schapter
}

\setlength{\parskip}{.5em} 					% Redefinição do espaço entre parágrafos

%==============================================================================
%%% Hack: Adicionando Bibliografia, Lista de Figuras, Lista de Tabelas,
%%% Lista de Símbolos e Lista de Abreviaturas ao Sumário
%==============================================================================
% Bibliografia
\renewenvironment{thebibliography}[1]{
	\chapter*{
		\bibname
		\refstepcounter{chapter}
		\addcontentsline{toc}{chapter}{\bibname}    % Thanks to D.Farfan
		\@mkboth{\MakeUppercase{\bibname}}
		{\MakeUppercase{\bibname}}
	}       										% Thanks to ffs :-)
	
	\singlespacing			
	\flushleft				
	\list{\@biblabel{\@arabic\c@enumiv}}%
	{
		\settowidth\labelwidth{\@biblabel{#1}}%
		\leftmargin\labelwidth
		\advance\leftmargin\labelsep
		\@openbib@code
		\usecounter{enumiv}%
		\let\p@enumiv\@empty
		\renewcommand\theenumiv{\@arabic\c@enumiv}
	}

	\sloppy
	\clubpenalty4000
	\@clubpenalty \clubpenalty
	\widowpenalty4000%Y
	\sfcode`\.\@m
}
{
	\def\@noitemerr
	{\@latex@warning{Empty `thebibliography' environment}}%
	\endlist
}

% Lista de Figuras
\renewcommand\listoffigures{
	%\begin{spacing}{1.5}								% Espaçamento padrão
	
	\if@twocolumn
		\@restonecoltrue\onecolumn
	\else
		\@restonecolfalse
	\fi
	
	\chapter*{
		\listfigurename
		\@mkboth{\MakeUppercase{\listfigurename}}
		{\MakeUppercase{\listfigurename}}
	}      											   % Thanks to Reinaldo Fortes

	%\addcontentsline{toc}{chapter}{\listfigurename}   % Thanks to D.Farfan
	\pdfbookmark[0]{\listfigurename}{listafiguras}
	\@starttoc{lof}
	\if@restonecol
		\twocolumn
	\fi
	
	%\end{spacing}										% Espaçamento padrão
}

% Lista de Tabelas
\renewcommand\listoftables{
	%\begin{spacing}{1.5}								% Espaçamento padrão
	
	\if@twocolumn
		\@restonecoltrue\onecolumn
	\else
		\@restonecolfalse
	\fi
	
	\chapter*{
		\listtablename
		\@mkboth{\MakeUppercase{\listtablename}}
		{\MakeUppercase{\listtablename}}
	}     												% Thanks to Reinaldo Fortes
	
	%    \addcontentsline{toc}{chapter}{\listtablename} % Thanks do D.Farfan
	
	\pdfbookmark[0]{\listtablename}{listatabelas}
	\@starttoc{lot}
	\if@restonecol\twocolumn\fi
	
	%\end{spacing}										% Espaçamento padrão
}

% Redefinição do Sumário
\renewcommand\tableofcontents{%
	%\begin{spacing}{1.5}								% Espaçamento padrão
	
	\if@twocolumn
		\@restonecoltrue\onecolumn
	\else
		\@restonecolfalse
	\fi
	
	\chapter*{
		\contentsname
		\pdfbookmark[0]{\namecontents}{contents}
		\@mkboth{\MakeUppercase{\contentsname}}
		{\MakeUppercase{\contentsname}}
	}       											% Thanks to Reinaldo Fortes
	
	\@starttoc{toc}
	\if@restonecol
		\twocolumn
	\fi
	
	%\end{spacing}										% Espaçamento padrão
}

% Lista de Símbolos
\newcommand\listofsymbols{
	%\begin{spacing}{1.5}								% Espaçamento padrão
	
	\if@twocolumn
		\@restonecoltrue\onecolumn
	\else
		\@restonecolfalse
	\fi
	
	\chapter*{
		\listofsymbolsname
		\@mkboth{\MakeUppercase{\listofsymbolsname}}
		{\MakeUppercase{\listofsymbolsname}}
	}													% Thanks to Reinaldo Fortes
	
	%    \addcontentsline{toc}{chapter}{\listofsymbolsname} % Thanks to ffs  :-)
	
	\pdfbookmark[0]{\listofsymbolsname}{listasimbolos}
	\@starttoc{lof}
	\if@restonecol
		\twocolumn
	\fi

	%\end{spacing}										% Espaçamento padrão
}

% Lista de Abreviaturas
\newcommand\listofabbreviations{
	%\begin{spacing}{1.5}								% Espaçamento padrão
	
	\begin{flushleft}
		\if@twocolumn
			\@restonecoltrue\onecolumn
		\else
			\@restonecolfalse
		\fi
		
		\chapter*{
			{\listofabbreviationsname}
			\@mkboth{\MakeUppercase{\listofabbreviationsname}}		% Thanks to Reinaldo Fortes
			{\MakeUppercase{\listofabbreviationsname}}
		}%
		
		%\addcontentsline{toc}{chapter}{\listofabbreviationsname}%  Thanks to ffs  :-)
		
		\pdfbookmark[0]{\listofabbreviationsname}{listaabreviaturas}
		\@starttoc{lof}%
		\if@restonecol
			\twocolumn
		\fi
	\end{flushleft}

	%\end{spacing}										% Espaçamento padrão
}


% Thanks to Dr. Alan Wilter Sousa da Silva - alan@lncc.br
\renewcommand{\addcontentsline}[3]{
	\begingroup
		\let\label\@gobble
		\let\textlatin\@firstofone
	
		\ifx\@currentHref\@empty
			\Hy@Warning{
				No destination for bookmark of \string\addcontentsline, \MessageBreak destination is added
			}
			\phantomsection
		\fi
		
		\expandafter\ifx\csname toclevel@#2\endcsname\relax
			\begingroup
				\def\Hy@tempa{#1}%
				
				\ifx\Hy@tempa\Hy@bookmarkstype
					\Hy@WarningNoLine{
						bookmark level for unknown #2 defaults to 0
					}
				\else
					\Hy@Info{
						bookmark level for unknown #2 defaults to 0
					}
				\fi
			\endgroup
			\expandafter\gdef\csname toclevel@#2\endcsname{0}%
		\fi
		
		\edef\Hy@toclevel{\csname toclevel@#2\endcsname}%
		\Hy@writebookmark{\csname the#2\endcsname}%
		{#3}%
		{\@currentHref}%
		{\Hy@toclevel}%
		{#1}%
		\ifHy@verbose
			\typeout{pdftex: bookmark at \the\inputlineno:
				{\csname the#2\endcsname}
				{#3}
				{\@currentHref}
				{\Hy@toclevel}
				{#1}
			}
		\fi
		
		\addtocontents{#1}{
			\protect\contentsline{#2}{#3}{\thepage}{\@currentHref}
		}
	\endgroup
}

%==============================================================================
%%% Cabeçalhos e Rodapés
%==============================================================================
\pagestyle{fancy}
\lhead{\leftmark}
\chead{}
\rhead{\thepage}
\lfoot{}\cfoot{}\rfoot{}
\renewcommand{\headrulewidth}{0.5pt}
\renewcommand{\footrulewidth}{0pt}

\if@draft								% Imprimir rodapé com versão preliminar
	\rfoot{\emph{\itadraftbabel\ \today}}
\fi

%==============================================================================
%%% Membros da banca examinadora
%==============================================================================
\newcounter{mecounter}
\newcommand{\examiner}[6]{
	\addtocounter{mecounter}{1}
	\expandafter\gdef\csname meCarrer\alph{mecounter}\endcsname{#1}
	\expandafter\gdef\csname meDegree\alph{mecounter}\endcsname{#2}
	\expandafter\gdef\csname meName\alph{mecounter}\endcsname{#3}
	\expandafter\gdef\csname meGender\alph{mecounter}\endcsname{#4}
	\expandafter\gdef\csname meTitle\alph{mecounter}\endcsname{\genderCarrer{#1}{#4}~\genderDegree{#2}{#4}}
	\expandafter\gdef\csname meRole\alph{mecounter}\endcsname{#5}
	\expandafter\gdef\csname meInst\alph{mecounter}\endcsname{#6}
}

\newcommand{\@memberslist}{
	\begin{tabular}{l l l l}
		\@ifundefined{meTitlea}{}{\meTitlea\@ & \meNamea\@ & \hspace{1cm} \meRolea\@ & - \hspace{0.5cm} \meInsta\@} \\
		
		\advisorTitle & \advisorName & \hspace{1cm} \advisorRole & - \hspace{0.5cm} \advisorInst \\
		
		\@ifundefined{@coadvisor}{}{\coadvisorTitle & \coadvisorName & \hspace{1cm} \coadvisorRole & - \hspace{0.5cm} \coadvisorInst \\}
		
		\@ifundefined{meTitleb}{}{\meTitleb\@ & \meNameb\@ & \hspace{1cm} \meRoleb\@ & - \hspace{0.5cm} \meInstb\@} \\
		
		\@ifundefined{meTitlec}{}{\meTitlec\@ & \meNamec\@ & \hspace{1cm} \meRolec\@ & - \hspace{0.5cm} \meInstc\@} \\
		
		\@ifundefined{meTitled}{}{\meTitled\@ & \meNamed\@ & \hspace{1cm} \meRoled\@ & - \hspace{0.5cm} \meInstd\@} \\
		
		\@ifundefined{meTitlee}{}{\meTitlee\@ & \meNamee\@ & \hspace{1cm} \meRolee\@ & - \hspace{0.5cm} \meInste\@} \\
		
		\@ifundefined{meTitlef}{}{\meTitlef\@ & \meNamef\@ & \hspace{1cm} \meRolef\@ & - \hspace{0.5cm} \meInstf\@} \\
		
		\@ifundefined{meTitleg}{}{\meTitleg\@ & \meNameg\@ & \hspace{1cm} \meRoleg\@ & - \hspace{0.5cm} \meInstg\@} \\
		
		\@ifundefined{meTitleh}{}{\meTitleh\@ & \meNameh\@ & \hspace{1cm} \meRoleh\@ & - \hspace{0.5cm} \meInsth\@}
	\end{tabular}
}

%==============================================================================
%%% Palavras-chave
%==============================================================================
\newcommand{\keywords}[1]{
	\vfill \noindent\textbf{Palavras--Chave:} #1
}

%==============================================================================
%%% Keywords
%==============================================================================
\newcommand{\englishkeywords}[1]{
	\vfill \noindent\textbf{Keywords:} #1
}

%==============================================================================
% Definição das palavras-chave para CIP (informado pela Biblioteca)
%==============================================================================
\newcounter{kwcipcounter}
\newcommand{\kwcip}[1]{
	\addtocounter{kwcipcounter}{1}
	\expandafter\gdef\csname kwcip\alph{kwcipcounter}\endcsname{#1}
}

\newcommand{\@cipkw}{
	\@ifundefined{kwcipa}{}{1.~\expandafter\MakeUppercase\kwcipa\@. }%
	\@ifundefined{kwcipb}{}{2.~\expandafter\MakeUppercase\kwcipb\@. }%
	\@ifundefined{kwcipc}{}{3.~\expandafter\MakeUppercase\kwcipc\@. }%
	\@ifundefined{kwcipd}{}{4.~\expandafter\MakeUppercase\kwcipd\@. }%
	\@ifundefined{kwcipe}{}{5.~\expandafter\MakeUppercase\kwcipe\@. }%
	\@ifundefined{kwcipf}{}{6.~\expandafter\MakeUppercase\kwcipf\@. }%
	\@ifundefined{kwcipg}{}{7.~\expandafter\MakeUppercase\kwcipg\@. }%
	\@ifundefined{kwciph}{}{8.~\expandafter\MakeUppercase\kwciph\@. }%
	\@ifundefined{kwcipi}{}{9.~\expandafter\MakeUppercase\kwcipi\@. }%
	\@ifundefined{kwcipj}{}{10.~\expandafter\MakeUppercase\kwcipj\@. }%
}

%==============================================================================
%%% Hack: Centralizando \chapter*
%%% \@makeschapterhead é chamado por \chapter*
%==============================================================================
\def\@makeschapterhead#1{
	\vspace*{50\p@}
	{
		\parindent \z@ \centering
		\normalfont
		\interlinepenalty\@M
		\Huge \bfseries #1\par\nobreak
		\vskip 40\p@
	}
}

%==============================================================================
%%% Redefinição do comando \appendix para mostrar 'Apêndice' ou 'Anexo' no Sumário
%==============================================================================
\renewcommand{\appendix}{
	
	\cleardoublepage
	\providecommand\anexoautorefname{appendix}
	\setcounter{chapter}{0}%
	\setcounter{section}{0}%
	\def\Hy@appendixstring{appendix}			% novo string no *.toc, para fazer o hyperlink funcionar
	\def\anexo{%
		\@ifundefined{chapter}%
		{\gdef\theHsection{\Alph{section}}}%
		{\gdef\theHchapter{\Alph{chapter}}}%
		
		\xdef\Hy@chapapp{\Hy@appendixstring}%
		\gdef\@chapapp{\appendixname}%
		
		\addtocontents{toc}{
			\SetAppendixTocName{\appendixname}
		}%
		
		\gdef\thechapter{\@Alph\c@chapter}%
	}
	\anexo
}

%==============================================================================
%%% Acrescenta o comando SetappendixTocName ao toc
%==============================================================================
\DeclareRobustCommand{\SetAppendixTocName}[1]{\gdef\AppendixTocName{#1}}

%==============================================================================
%%% Novo comando: \annex
%==============================================================================
\newcommand{\annex}{%
  \cleardoublepage
  \providecommand\anexoautorefname{annex}
	\setcounter{chapter}{0}%
	\setcounter{section}{0}%
	\def\Hy@appendixstring{annex}% %%novo string no *.toc, para fazer o hyperlink funcionar
	\def\anexo{%
  \@ifundefined{chapter}%
    {\gdef\theHsection{\alph{section}}}%
    {\gdef\theHchapter{\alph{chapter}}}%
  \xdef\Hy@chapapp{\Hy@appendixstring}%

  \gdef\@chapapp{\nameannex}%
  \addtocontents{toc}{\SetAppendixTocName{\nameannex}}%
  \gdef\thechapter{\@Alph\c@chapter}%
  }
  \anexo
}


%==============================================================================
%%% Redefinição de mainmatter para numeração continua de páginas
%==============================================================================
\renewcommand\mainmatter{
	\cleardoublepage
	\@mainmattertrue
	\renewcommand{\thepage}{\arabic{page}}
}

%==============================================================================
%%% 'Desligar' hifenação...
%==============================================================================
\newlanguage\nohyphen
\newcommand\nh[1]{{\language\nohyphen #1}}

%==============================================================================
%%% Redefinição de l@chapter, l@section e l@subsection para acrescentar
%%% a palavra 'Anexo' ou 'Apêndice' no toc e ajustar identação para
%%% "nada identado" ... :-)
%==============================================================================
\renewcommand*\l@chapter[2]{
	\ifnum \c@tocdepth >\m@ne
	\addpenalty{-\@highpenalty}%
	\vskip -3ex \@plus\p@
	\setlength\@tempdima{1.5em}%
	\begingroup
		\ifx\AppendixTocName\undefined
		\else
			\renewcommand{\numberline}[1]{\AppendixTocName~##1~~--~~}%
		\fi
		
		\parindent \z@ \rightskip \@pnumwidth
		\parfillskip -\@pnumwidth
		\leavevmode \bfseries
		\advance\leftskip\@tempdima
		\hskip -\leftskip
		\@dottedtocline{0}{0em}{1.8em}{\normalfont\scshape\large{#1}}{#2}
		\penalty\@highpenalty
	\endgroup
	\fi
}

% l@section
\renewcommand*{\l@section}[2]{
	\@dottedtocline{1}{0.5cm}{2.8em}{\textbf{#1}}{\hss #2}
}

% l@subsection
\renewcommand*{\l@subsection}{
	\@dottedtocline{2}{1cm}{3.7em}
}

%==============================================================================
%%% Numerar o subsubsection conforme definido anteriormente
%==============================================================================
\setcounter{secnumdepth}{4}

%==============================================================================
%%% Glossário
%==============================================================================
\newcommand{\itaglossary}{
	\cleardoublepage
	\thispagestyle{empty}%
	
	\@mkboth{\MakeUppercase{\itaglossarynamebabel}}%
	        {\MakeUppercase{\itaglossarynamebabel}}%
	
	\gdef\glossaryname{\itaglossarynamebabel}
	\refstepcounter{chapter}
}

%==============================================================================
% %%% Centralizar e ajustar os espaçamentos dos captions das figuras e tabelas 
%==============================================================================
\def\fnum@figure{\centering\MakeUppercase{\figurename}\nobreakspace\thefigure}
\def\fnum@table{\centering\MakeUppercase{\tablename}\nobreakspace\thetable}
\setlength\abovecaptionskip{6\p@}
\setlength\belowcaptionskip{4\p@}

% Redefinição do caracter que separa o label Figura ou Tabela do texto
%\renewcommand*\captionlabeldelim{\ --}

%==============================================================================
% %%% Indentar a segunda linha de captions longos das listas 
%     de figuras e tabelas
%==============================================================================
\renewcommand*{\l@figure}{
	\@dottedtocline{1}{0cm}{1cm}
}

\def\FigureITA@dottedtocline#1#2#3#4#5{
	\begingroup
		\def\numberline##1{\hb@xt@\@tempdima{\itafigurenamebabel~##1 --\hfil}}%
  		\@dottedtocline{#1}{#2}{#3}{#4}{#5}%
	\endgroup
}

\def\TableITA@dottedtocline#1#2#3#4#5{
	\begingroup
	\def\numberline##1{\hb@xt@\@tempdima{\itatablenamebabel~##1 --\hfil}}%
	\@dottedtocline{#1}{#2}{#3}{#4}{#5}%
	\endgroup
}

% faz a lista de figuras e de tabelas usarem as novas versões dos \dottedtocline
\renewcommand*\l@figure{\FigureITA@dottedtocline{1}{0cm}{2.8cm}}
\renewcommand*\l@table{\TableITA@dottedtocline{1}{0cm}{2.8cm}}


%==============================================================================
% Folha de Registro do Documento
%==============================================================================
\newcommand{\esp}{.2cm}

% Definições dos valores dos campos do formulário
\newcommand{\FRDitadocnro}[1]{\gdef\FRDitadocnrobabel{#1}}
\newcommand{\FRDitaorgaointerno}{\FRDitaorgaointernobabel}
\newcommand{\FRDitapalavrasautor}[1]{\gdef\FRDitapalavrasautorbabel{#1}}
\newcommand{\FRDitapalavrasresult}[1]{\gdef\FRDitapalavrasresultbabel{#1}}
\newcommand{\FRDitapalavraapresentacao}{\FRDitapalavraapresentacaobabel}
\newcommand{\FRDitaresumo}[1]{\gdef\FDRitaabstractcontent{#1}}

%  Opções para: 
%  Nacional ou Internacional -- N/I
%  Ostensivo, Reservado, Confidencial ou Secreto -- O/R/S
\newcommand{\FRDitaOpcoes}[2]{
	\gdef\N{\ } \gdef\I{\ }      
	\gdef\O{\ } \gdef\R{\ } \gdef\C{\ } \gdef\S{\ }
	
	\ifthenelse{\equal{#1}{N}}{\gdef\N{\textbf{X}}}{}
	\ifthenelse{\equal{#1}{I}}{\gdef\I{\textbf{X}}}{}
	
	\ifthenelse{\equal{#2}{O}}{\gdef\O{\textbf{X}}}{}
	\ifthenelse{\equal{#2}{R}}{\gdef\R{\textbf{X}}}{}
	\ifthenelse{\equal{#2}{S}}{\gdef\S{\textbf{X}}}{}
}

% Folha de Registro de Documento
\newcommand{\itaFRD}{
	
	\if@openright
		\cleardoublepage
	\else
		\clearpage
	\fi

	% Redefinição das margens
	\setlength{\hoffset}{-1in}        % remove padrão do LaTeX
	\setlength{\voffset}{-1in}        % remove padrão do LaTeX
	\setlength{\oddsidemargin}{1cm}   % margem esquerda, páginas ímpares
	\setlength{\evensidemargin}{1cm}  % margem esquerda, páginas pares
	\setlength{\topmargin}{3cm}       % margem superior
	\newlength{\bottommargin}         % margem inferior é deduzida em LaTeX
	\setlength{\bottommargin}{2cm}    % margem inferior
	\setlength{\headheight}{0pt}      % altura do cabeçalho
	\setlength{\headsep}{0pt}         % espaço entre cabeçalho e texto
	\setlength{\footskip}{0mm}        % espaço entre rodapé e texto
	\setlength\textheight{            % altura disponível para texto
		\paperheight-\topmargin-\bottommargin-\headheight-\headsep-\footskip
	}
	
	% margem direita é deduzida em LaTeX... alterando isso...
	\newlength{\rightmarginn}        
	\setlength{\rightmarginn}{2cm}    % margem direita
	\setlength\textwidth{             % largura disponível para texto
		\paperwidth-\oddsidemargin-\rightmarginn
	}

	\singlespacing
	\footnotesize
	\setlength{\baselineskip}{0.1\baselineskip}
	\setlength{\parskip}{.1em}
	\setlength{\extrarowheight}{2pt}

	\lhead{}\rhead{}
	
	\thispagestyle{empty}
	\pdfbookmark[0]{\itapagedocregbabel}{bla}
	
	\begin{tabularx}{\textwidth}{|X|X|X|X|}
		\hline
		\multicolumn{ 4}{|c|}{} \\
		\multicolumn{ 4}{|c|}{\normalsize FOLHA DE REGISTRO DO DOCUMENTO} \\
		\multicolumn{ 4}{|c|}{} \\
		\hline
		
		$^{1.}$ \scriptsize CLASSIFICAÇÃO/TIPO: & 
		$^{2.}$ \scriptsize DATA:               & 
		$^{3.}$ \scriptsize DOCUMENTO Nº:       & 
		$^{4.}$ \scriptsize Nº DE PÁGINAS:             \tabularnewline 
		
		\centering{
			\ifthenelse{\equal{\itaoptionlevel}{tg}}
			{TC}
			{\ifthenelse{\equal{\itaoptionlevel}{msc}}
			{DM}
			{TD}
			}
		}  & 
		
		\centering{\DTMusedate{publDate}}   & 
		
		\centering\scriptsize \FRDitadocnrobabel  &  
		
		\centering\NumPenultimaPagina                  \tabularnewline[\esp]
		\hline
		
		\multicolumn{ 4}{|p{17cm}|}{$^{5.}$ \scriptsize TÍTULO E SUBTÍTULO:} \\
		\multicolumn{ 4}{|l|}{\nh{\parbox[b]{.95\textwidth}{\@title}}} \\[\esp]
		\hline
		\multicolumn{ 4}{|l|}{$^{6.}$ \scriptsize AUTOR\ifthenelse{\equal{\@authorGender}{masc}}{}{A}(ES):} \\
		\multicolumn{ 4}{|l|}{\textbf{\@authorFirst\ \@authorLast}} \\[\esp]
		\hline
		\multicolumn{ 4}{|l|}{$^{7.}$ \scriptsize INSTITUIÇÃO(ÕES)/ÓRGÃO(S) 
		INTERNO(S)/DIVISÃO(ÕES):} \\
		\multicolumn{ 4}{|l|}{\FRDitaorgaointernobabel} \\[\esp]
		\hline
		\multicolumn{ 4}{|l|}{$^{8.}$ \scriptsize PALAVRAS-CHAVE SUGERIDAS PEL\ifthenelse{\equal{\@authorGender}{masc}}{O}{A} AUTOR\ifthenelse{\equal{\@authorGender}{masc}}{}{A}:} \\
		\multicolumn{ 4}{|p{17cm}|}{\FRDitapalavrasautorbabel} \\[\esp]
		\hline
		\multicolumn{ 4}{|l|}{$^{9.}$ \scriptsize PALAVRAS-CHAVE RESULTANTES DE INDEXAÇÃO:} \\
		\multicolumn{ 4}{|p{17cm}|}{\FRDitapalavrasresultbabel} \\[\esp]
		\hline
		\multicolumn{ 4}{|l|}{$^{10.}$ {\scriptsize APRESENTAÇÃO:}  \hfill 
		(\N) \textbf{Nacional} \ \ \ 
		(\I) \textbf{Internacional} } \\
		\multicolumn{ 4}{|p{17cm}|}{\DTMsetdatestyle{ddmmyyyy} \DTMsetup{datesep=/} \ifthenelse{\equal{\itaoptionlevel}{tg}}{\FRDitapalavraapresentacaotgbabel}{\FRDitapalavraapresentacaobabel}} \\[\esp]
		\hline
		\multicolumn{ 4}{|l|}{$^{11.}$ \scriptsize RESUMO:} \\
		\multicolumn{ 4}{|c|}{\parbox[t][12cm]{.95\textwidth}{\FDRitaabstractcontent}} \\[\esp]
		\hline
		\multicolumn{ 4}{|l|}{$^{12.}$ \scriptsize GRAU DE SIGILO:} \\
		\multicolumn{ 4}{| l|}{\hspace{2.0cm} (\O) \textbf{OSTENSIVO} \hspace{2.0cm} 
		(\R) \textbf{RESERVADO} \hspace{2.0cm} 
		(\S) \textbf{SECRETO} \hspace{2.0cm}} \\[\esp]
		\hline
     \end{tabularx}
}


%==============================================================================
%%% Hack: Retirando a palavra Capitulo mas mantendo Anexo e Apendice...
%%% \@makechapterhead é chamado por \chapter
%==============================================================================
\def\@makechapterhead#1{
	\vspace*{50\p@}
	
	{
		\parindent \z@ \raggedright \normalfont
		\ifnum \c@secnumdepth >\m@ne
		
		\if@mainmatter
			% se "anexo" nao foi definido, altera cabecalho do capitulo..
			% senao, mantem a palavra Anexo ou Apendice
			\@ifundefined{anexo}      
				{\huge\bfseries \thechapter\space}
				{\huge\bfseries \@chapapp\space \thechapter \space -}
			\fi
		\fi
		
		\interlinepenalty\@M
		\Huge \bfseries #1\par\nobreak
		\vskip 40\p@
	}
}
  
%==============================================================================
%%% Obtendo o número da penúltima página, visto que a FRD nao é contada ...
%%% Thanks to Dr. Miguel V. S. Frasson!! :-)
%==============================================================================
\def\deixasosegundo#1#2#3!{#2}
\makeatletter
\AtBeginDocument{
	% se a referência a última página não existe, definir como ??
	\expandafter\ifx
		\csname r@LastPage\endcsname\relax
		\newcommand\NumPenultimaPagina{{\bfseries??}}%
	\else
		\newcounter{butlastpage}%
		\setcounter{butlastpage}{\expandafter\deixasosegundo\r@LastPage!}%
		\addtocounter{butlastpage}{-1}%
		\newcommand\NumPenultimaPagina{}%
		\edef\NumPenultimaPagina{\arabic{butlastpage}}%
	\fi
}

\endinput