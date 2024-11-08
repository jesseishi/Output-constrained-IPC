# IEA Wind 15-MW RWT: OpenFAST Model

Adapted from
[GitHub: IEA 15
MW](https://github.com/IEAWindTask37/IEA-15-240-RWT/tree/master/openFAST).
Changes:
- Remove unused files when using the monopile configuration and remove 'Monopile'
  specification from input file names.
- Disabled all echo's and summaries.
- Set control modes to 4 in ServoDyn so we can control the model from Simulink.
- Set `NumCrctn` to 1.
- Set `OutFileFmt` to 2 to save data by only saving the binary file and not the text
  file too.
- Remove some outputs from the different OutLists to save data.
- Add some wind files.
