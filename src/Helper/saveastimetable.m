function saveastimetable(fname, simOut, OutList, fast)
% Compile the logged signals. We have two sources for this: 
% logging from Simulinks and the outb file from openFAST.
TT1 = extractTimetable(simOut.logsout);
% With the default format (x sec), readtimetable doesn't work. So
% change the format so we can read this timetable gain. This also 
% makes it friendly to load with other programs such as Python.
TT1.Time.Format = 'hh:mm:ss.SSS';

TT2 = out2TT(simOut.OutData, OutList);
TT = [TT1, TT2];

% Save the logged signals.
writetimetable(TT, fname);

% We should now have a `.outb` file in the same location as the
% FAST_InputFile. However, we can use multiple controller with the same
% input file, so we also move it to the outDir.
% I keep this in outb format so that it is also readable using e.g. Python.
outbFileName = replace(fast.FAST_InputFile, '.fst', '.SFunc.outb');
outbFile = fullfile(fast.FAST_directory, outbFileName);
outbFileNew = replace(fname, '.csv', '.SFunc.outb');
movefile(outbFile, outbFileNew)
end