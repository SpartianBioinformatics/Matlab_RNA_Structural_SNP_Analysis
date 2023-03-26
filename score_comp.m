function [ score ] = score_comp( str1, str2 )
%SCORE_COMP compares two RNA structures in dot-bracket notation
%   Takes two dot-bracket notation RNA secondary structures as input and
%   compares their structure by aligning them and outputting a score as a
%   function of the RNA's length. The highest score of 1 indicates the same
%   structure and lower scores (lowest is 0) indicate dissimilar
%   structures.

str1 = strrep(str1,'.','A'); str1 = strrep(str1,'(','C'); str1 = strrep(str1,')','G');
str2 = strrep(str2,'.','A'); str2 = strrep(str2,'(','C'); str2 = strrep(str2,')','G');

score = nwalign(str1,str2,'Alphabet','NT','ScoringMatrix',eye(4),'GAPOPEN',.5,'EXTENDGAP',0.25);
score = score / max(length(str1),length(str2));

end
