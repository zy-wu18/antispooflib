classdef Logger < handle
    % LOGGER: To facilitate logging management
    %   IndStack: indent stack, non-negative integer
    %   IndStyle: indent style, default = '  '
    
    properties
        Enabled
        IndStack
        IndStyle
        BarLength
        BarStyle
    end

    properties(Hidden)
        Bar
        BarFullLog
        BarProg
        Tab
        TabFullLog
        TabColumnNum
        TabLineNum
        TabLineFormat
        TabLineName
    end
    
    methods
        %% Construct and initialize an instance of Logger
        function obj = Logger()
            obj.Enabled = true;
            obj.IndStack = 0;
            obj.IndStyle = '  ';
            obj.BarLength = 15;
            obj.BarStyle = '>.';
            obj.BarProg = 0;
            obj.Bar = repmat(obj.BarStyle(2), [1, obj.BarLength]);
            obj.BarFullLog = [];
            obj.Tab = [];
            obj.TabColumnNum = 0;
            obj.TabLineNum = 0;
            obj.TabFullLog = [];
            obj.TabLineName = {};
            obj.TabLineFormat = {};
        end
        
        %% Basic logging
        function writeLine(obj, varargin)
            narginchk(1, inf);
            varargin{1} = [char(varargin{1}), '\n'];
            obj.write(varargin{1:end});
        end
        
        function write(obj, varargin)
            if(~obj.Enabled)
                return;
            end
            narginchk(1, inf);
            fprintf(repmat(obj.IndStyle, [1, obj.IndStack]));
            if(length(varargin) < 2)
                fprintf(char(varargin{1}));
            else
                fprintf(char(varargin{1}), varargin{2:end}); 
            end
        end

        %% Operations on indent stack
        function enStack(obj, varargin)
            if(nargin == 0)
                varargin = cell([]);
            end
            if(~isempty(varargin))
                obj.writeLine(varargin{1:end});
            end
            obj.IndStack = obj.IndStack + 1;
        end

        function deStack(obj, varargin)
            obj.IndStack = obj.IndStack - 1;
            if(nargin == 0)
                varargin = cell([]);
            end
            if(~isempty(varargin))
                obj.writeLine(varargin{1:end});
            end
        end
        
        %% Progression bar logging
        function refreshBar(obj, val_cur, val_all, bar_format)
            if(nargin == 3)
                bar_format = 'Progress: %s %5d/%5d';
            end
            if(obj.BarLength*val_cur/val_all > obj.BarProg)
                step = max(1, floor(obj.BarLength/val_all));
                obj.Bar = [repmat(obj.BarStyle(1), [1, step]), obj.Bar(1:end-step)];
                if(obj.BarProg > 0)
                    fprintf(repmat('\b', [1, length(obj.IndStyle)*obj.IndStack]));
                end
                fprintf(repmat('\b', [1, length(obj.BarFullLog)]));
                obj.BarFullLog = sprintf(bar_format, obj.Bar, val_cur, val_all);
                obj.write(obj.BarFullLog);
                obj.BarProg = obj.BarProg + step;
            end
        end

        function resetBar(obj)
            fprintf('\n');
            obj.BarProg = 0;
            obj.Bar = repmat(obj.BarStyle(2), [1, obj.BarLength]);
            obj.BarFullLog = [];
        end

        %% Tablet logging
        function setTable(obj, lformat, lname)
            obj.resetTable();
            if(length(lname)~=length(lformat))
                obj.write("[Error] Logger.setTableFormat requires line_format and line_name with the same length, ");
                obj.writeLine("but received cells with length %d and %d.", length(lname), length(lformat));
                return;
            end
            obj.TabLineNum = length(lformat);
            obj.TabLineFormat = string(lformat);
            lname_maxlen = 0;
            for l = 1:obj.TabLineNum
                if(length(lname{l}) > lname_maxlen)
                    lname_maxlen = length(lname{l});
                end
            end
            for l = 1:obj.TabLineNum
                obj.TabLineName{l} = [lname{l}, repmat(' ', [1, lname_maxlen-length(lname{l})])];
            end
        end

        function refreshTable(obj, new_col)
            if(size(new_col, 1) ~= obj.TabLineNum)
                obj.write("[Error] Logger.refreshTable requires new column with size (%d, ~), ", obj.TabLineNum);
                obj.writeLine("but received one with size (%d, %d)", size(new_col, 1), size(new_col, 2));
                return;
            end
            obj.Tab = [obj.Tab, new_col];
            if(obj.TabColumnNum > 0)
                fprintf(repmat('\b', [1, length(obj.IndStyle)*obj.IndStack*obj.TabLineNum]));
                for l = 1:obj.TabLineNum
                    fprintf(repmat('\b', [1, length(obj.TabFullLog{l})+1]));
                end
            end
            tab_full_log = repmat("", [1, obj.TabLineNum]);
            for l = 1:obj.TabLineNum
                tab_full_log{l} = char(strcat(string(obj.TabLineName{l}), " = "));
                for c = 1:obj.TabColumnNum + size(new_col, 2)
                    tab_full_log_tmp = strcat(tab_full_log{l}, sprintf(string(obj.TabLineFormat{l}), obj.Tab(l, c)));
                    tab_full_log{l} = char(tab_full_log_tmp);
                end
                obj.writeLine(tab_full_log{l});
            end

            obj.TabFullLog = tab_full_log;
            obj.TabColumnNum = obj.TabColumnNum + size(new_col, 2);
        end

        function resetTable(obj)
            obj.Tab = [];
            obj.TabColumnNum = 0;
            obj.TabLineNum = 0;
            obj.TabFullLog = [];
            obj.TabLineName = {};
            obj.TabLineFormat = {};
        end
    end
end
