%%%%%
% Mathias' LaTeX class for typesetting notes
%
% Copy this class to $TEXMF/tex/latex/notes/notes.cls
% and do "\documentclass{notes}" in the pre-amble.
% 
% To find $TEXMF (which is usuall ~/texmf), do
%
%   kpsewhich -var-value=TEXMFHOME
%
% To check if LaTeX will find the class when
% compiling, run
%
%   kpsewhich notes.cls
%
% which should show ~/texmf/tex/latex/notes/notes.cls
%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{notes}[2017/05/24 Notes]

\LoadClass[12pt,a4paper]{article}

\RequirePackage[utf8]{inputenc}
\RequirePackage[T1]{fontenc}
\RequirePackage{geometry}
\RequirePackage{fancyhdr,lastpage}
\RequirePackage{amsmath}
\RequirePackage{graphicx}
\RequirePackage[labelfont=bf]{caption}
\RequirePackage{color}
\RequirePackage{comment}
\RequirePackage{titlesec}

% Allow to use fancier font if available
%\AtBeginDocument{
\IfFileExists{tgbonum.sty}{
    \RequirePackage{tgbonum}
    %\fontfamily{cmr}\selectfont
    \renewcommand{\rmdefault}{lmr}
    %\renewcommand{\familydefault}{\rmdefault}
}{}
%}

\DeclareOption{altfont}{
    %\fontfamily{qbk}\selectfont
    \renewcommand{\rmdefault}{qbk}
    \renewcommand{\sfdefault}{lmss}

    \renewcommand*\maketitle{
        \begingroup
        \centering
        \sffamily\@maketitle
        \par\endgroup
    }
    \titleformat{\section}
        {\normalfont\sffamily\Large}
        {\thesection}{1em}{}
}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

%\ExecuteOptions{altfont}
\ProcessOptions\relax

\geometry{left=25mm,right=25mm,bindingoffset=0mm,top=20mm,bottom=20mm}

\makeatletter
\def\@maketitle{
	\newpage
	\null
	\noindent
	{\bfseries\LARGE\@title}\\
	{\small by \@author\hfill\@date}
    \vspace{2mm}
	\hrule\hfill
}
\makeatother

\newcommand{\DREAM}{\textsc{Dream}}

\newcommand{\bb}[1]{\boldsymbol{#1}}
\newcommand{\dd}{\mathrm{d}}
\newcommand{\thetap}{\theta_\mathrm{p}}

\newcommand{\xhat}{\hat{\bb{x}}}
\newcommand{\yhat}{\hat{\bb{y}}}
\newcommand{\zhat}{\hat{\bb{z}}}

\newcommand{\red}[1]{\textbf{\textcolor{red}{#1}}}

