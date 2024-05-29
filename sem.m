function [output] = sem (data) 
output=nanstd(data)/sqrt(size(data,1)); %Standard error of the mean
end