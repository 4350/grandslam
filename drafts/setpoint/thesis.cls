%%
%% Master's Thesis Class File
%%
%% Adapted from the PracTeX Journal class to follow SSE thesis guidelines
%% and in general be more appropriate. See original here:
%% http://tug.org/pracjourn/styles/latex/pracjourn.cls
%%
%%     This program is free software; you can redistribute it and/or modify
%%     it under the terms of the GNU General Public License as published by
%%     the Free Software Foundation; either version 2 of the License, or
%%     (at your option) any later version.
%% 
%%     This program is distributed in the hope that it will be useful,
%%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%     GNU General Public License for more details.
%% 
%%     You should have received a copy of the GNU General Public License
%%     along with this program; if not, write to the Free Software
%%     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%% 
\NeedsTeXFormat{LaTeX2e}[1999/12/01]
\ProvidesClass{thesis}
    [2016/11/19 v0.1n SSE MSc Thesis Class]

%% BASIC PAPER OPTIONS
% Note: change `draft' to `final' to include figures
\let\ExecuteOptions@ltx\ExecuteOptions
\def\ExecuteOptions#1{%
 \ExecuteOptions@ltx{a4paper,12pt,oneside,onecolumn,final}%
}%
\LoadClassWithOptions{article}
\let\ExecuteOptions\ExecuteOptions@ltx

\RequirePackage[a4paper,margin=2.50cm,footskip=1.25cm]{geometry}

%% GENERAL PACKAGES AND OPTIONS
% Colors and graphics; makes all links blue
\RequirePackage{color}
\RequirePackage{hyperref}
\RequirePackage{url}
\RequirePackage{graphicx}
\definecolor{linkcolour}{rgb}{0,0,0}
\hypersetup{colorlinks,breaklinks,
            linkcolor=linkcolour,citecolor=linkcolour,
            filecolor=linkcolour,urlcolor=linkcolour}

% We are Americans
\RequirePackage{csquotes}
\RequirePackage[american]{babel}

\RequirePackage{amsmath}
\RequirePackage{amssymb}

%% BIBLIOGRAPHY OPTIONS
% BIBLIOGRAPHY
\RequirePackage[
  hyperref=true,    % Make references hyperlinks
  style=authoryear  % Select style
]{biblatex}

% Put bibliography on its own page
\defbibenvironment{bibliography}
  {\begin{singlespacing}\list{}
    {\setlength{\leftmargin}{\bibhang}%
    \setlength{\itemindent}{-\leftmargin}%
    \setlength{\itemsep}{\bibitemsep}%
    \setlength{\parsep}{\bibparsep}}}
  {\endlist\end{singlespacing}}
  {\item}

% Make the bibliography heading an un-numbered section (not included in
% toc); refname uses the right translation, probably.
\defbibheading{bibliography}[\refname]{\section*{#1}}

\setlength\bibhang{0pt}               % Don't indent bibliography second line
\setlength\bibitemsep{\baselineskip}  % Separate entries by one line

%% TYPOGRAPHY OPTIONS

\RequirePackage{microtype}  % Black magic as far as I'm concerned
\RequirePackage{mathspec}   % Change font (only with XeTeX)
\RequirePackage{titlesec}   % Set section titles easily
\RequirePackage{titleps}    % For setting page headers/footers
\RequirePackage{enumerate}  % For enumerate environment

% Configures Minion Pro for everything with 1.5 line height
\setmainfont{Minion Pro}
\setmathsfont(Digits,Greek,Latin)[Numbers={Proportional}]{Minion Pro}
\setmathrm{Minion Pro}
\RequirePackage{setspace}
\onehalfspacing

%% SECTIONING
\setcounter{tocdepth}{2}         % Include only to subsection level
\numberwithin{equation}{section} % Number equations by sections

% Add "Appendix" section to ToC and start of appendices
\newcommand*{\Appendixautorefname}{Appendix}
\RequirePackage[toc,title]{appendix} % Ensure "Appendices" is added to ToC

%% TABLES AND FIGURES
% This puts table captions above tables and makes the "Table 1" label bold
\RequirePackage{caption}
\RequirePackage{subcaption}    % Multipage figures
\RequirePackage{booktabs}      % Nice tables
\RequirePackage{tabularx}      % Full-width tables
\RequirePackage{multirow}      % Span multiple rows
\RequirePackage{dcolumn}       % Aligned columns as below
\newcolumntype{d}{D{.}{.}{-1}} % Use `d` column for numerical columns
\newcolumntype{i}{D{.}{.}{0}}

\RequirePackage{placeins}      % \FloatBar­rier to stop floats
\RequirePackage{float}         % Better placement of floats

\captionsetup{font={bf,normalsize},skip=15pt}
\captionsetup[table]{position=top}
\captionsetup[figure]{belowskip=15pt,position=below}

% Environment for table descriptions
\newenvironment{longcaption}{\begin{minipage}{\textwidth}}{\\\end{minipage}}
