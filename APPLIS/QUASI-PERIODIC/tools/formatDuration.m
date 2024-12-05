function st = formatDuration(d,format)
% st = formatDuration(d,format)
% format is either 'standard' (default) or 'iso8601'.

if nargin<2
    format = 'standard' ;
end

if strcmp(format,'iso8601')
    sep = {'P','T','H','M','.','S'} ;
elseif strcmp(format,'standard')
    sep = {'','',':',':','.','s'} ;
else
    error('Unknown format requested')
end

stDays = '' ;
stHours = '' ;
stMinutes = '' ;

stCentiseconds = sprintf('%.2i',round(100*rem(d,1))) ;

stSeconds = sprintf('%.2i',mod(floor(d),60)) ;

if d>=60
    stMinutes = sprintf('%.2i',mod(fix(d/60),60)) ;
else
    sep{4} = '' ;
end

if d>=3600
    if strcmp(format,'iso8601')
        stHours = sprintf('%.2i',mod(fix(d/3600),24)) ;
    else
        stHours = sprintf('%.2i',fix(d/3600)) ;
    end
else
    sep{3} = '' ;
end

if d>=86400 && strcmp(format,'iso8601')
    stDays = sprintf('%iD',fix(d/86400)) ;
else
    sep{2} = '' ;
end

st = [sep{1},stDays,sep{2},stHours,sep{3},stMinutes,sep{4},stSeconds,...
    sep{5},stCentiseconds,sep{6}] ;
end

% Legacy : previous code
% switch format
%     case 'standard'
%         
%         % seconds
%         stSeconds = sprintf('%.2i',mod(floor(d),60)) ;
%         
%         % minutes
%         if d>=60
%             stMinutes = sprintf('%.2i:',mod(fix(d/60),60)) ;
%         else
%             stMinutes = '' ;
%         end
%         
%         % hours
%         if d>=3600
%             stHours = sprintf('%.2i:',fix(d/3600)) ;
%         else
%             stHours = '' ;
%         end
%         
%         % Final
%         st = [stHours,stMinutes,stSeconds] ;
%     
%     case 'iso8601' % 'precise' format : no year nor month
%         
%         % centiseconds
%         centiseconds = round(100*rem(d,1)) ;
%         if centiseconds>0
%             stCentiseconds = sprintf('.%.2iS',centiseconds) ;
%         else
%             stCentiseconds = 'S' ;
%         end
%         stSeconds = sprintf('%.2i',mod(floor(d),60)) ;
%         if d>=60
%             stMinutes = sprintf('%.2iM',mod(fix(d/60),60)) ;
%         else
%             stMinutes = 'M' ;
%         end
%         if d>=3600
%             stHours = sprintf('%.2iH',fix(d/3600)) ;
%         else
%             stHours = '' ;
%         end
%         if d>=86400
%             stDays = sprintf('%iDT',fix(d/86400)) ;
%         else
%             stDays = 'T' ;
%         end
%         st = ['P',stDays,stHours,stMinutes,stSeconds,stCentiseconds] ;
%     otherwise
%         error('formatDuration : unknown format')
% end