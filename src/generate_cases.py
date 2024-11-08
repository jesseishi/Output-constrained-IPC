# Make different cases to be studied with our constrained IPC controller.
# Imports.
import itertools
import os
import subprocess

# I like pathlib but the openfast toolbox uses os.path so it doesn't always combine very nicely.
from pathlib import Path
from typing import Any, Iterable, Optional

import numpy as np
import openfast_toolbox.case_generation.case_gen as case_gen
from openfast_toolbox.io.fast_input_file import FASTInputFile


def radps2rpm(radps: float) -> float:
    return radps * 60 / (2 * np.pi)


# TODO: I'm not super happy with this yet. The problem is that you want to use TurbSim
# to turn .dat files into .bts files and copy inputFiles using templateReplace. If you
# first run TurbSim you cannot let templateReplace delete unused files (which would be
# nice), but then you of course need to store your .dat files in the templateDir rather
# than the outputDir so you get quite a lot of .dat files in your templateDir. Nja,
# that's also not the worst thing ever.
def generate_cases(
    windfile: str,
    DOF: str = "full",
    randomSeeds: Optional[Iterable] = None,
    turbulenceIntencities: Optional[Iterable] = None,
    URefs: Optional[Iterable] = None,
    PLExps: Optional[Iterable] = None,
    TMaxs: Optional[Iterable] = None,
) -> None:
    # Set some defaults.
    randomSeeds = randomSeeds or [0]
    turbulenceIntencities = turbulenceIntencities or [0]
    URefs = URefs or [15]
    PLExps = PLExps or [0.1]
    TMaxs = TMaxs or [700]

    # Define some parameters.
    templateDir = Path("IEA-15-240-RWT")
    main_file = templateDir / "IEA-15-240-RWT.fst"
    outputDir = Path("Results/InputFiles")  # Will be created.

    # Base changes to the template .fst file.
    base_dict: dict[str, Any] = dict()

    # I'm not so interested in hydrodynamic loads for now.
    base_dict["CompHydro"] = 0

    # Enable/disable some DOFs.
    match DOF:
        case "Flap12":
            base_dict = case_gen.paramsStiff(base_dict)
            base_dict["EDFile|FlapDOF1"] = "True"
            base_dict["EDFile|FlapDOF2"] = "True"
            base_dict["EDFile|GenDOF"] = "True"
        case "FullDOF":
            pass
        case _:
            raise Exception(f"Don't know DOF {DOF}.")

    # This is not very modular but is ok for now, here I set some initial conditions
    # based on similar simulations.
    base_dict["EDFile|OoPDefl"] = 5.0
    base_dict["EDFile|IPDefl"] = -1.0
    base_dict["EDFile|BlPitch(1)"] = 11.6
    base_dict["EDFile|BlPitch(2)"] = 11.6
    base_dict["EDFile|BlPitch(3)"] = 11.6
    base_dict["EDFile|RotSpeed"] = 7.56
    base_dict["EDFile|TTDspFA"] = 0.17
    base_dict["EDFile|TTDspSS"] = -0.125

    # Collect all different parameter instances in a list.
    params = []
    for (
        randomSeed,
        turbulenceIntencity,
        URef,
        PLExp,
        TMax,
    ) in itertools.product(
        randomSeeds,
        turbulenceIntencities,
        URefs,
        PLExps,
        TMaxs,
    ):
        case_dict = base_dict.copy()

        # Process the wind file.
        windfilename, extension = windfile.split(".")
        # caseName = makeCaseName(DOF, randomSeed, turbulenceIntencity, URef, PLExp)
        caseName = makeName(
            "",
            DOF,
            "U",
            URef,
            "TI",
            turbulenceIntencity,
            "PLExp",
            PLExp,
            "seed",
            randomSeed,
            "T",
            TMax,
        )

        match extension:
            case "wnd":
                case_dict["InflowFile|WindType"] = 2
                case_dict["InflowFile|Filename_Uni"] = f"Wind\\{windfile}"
                # TODO: Overwrite the caseName. I made it more with the bts file in
                # mind. I should refactor this to be more general so that I can also use
                # steady wind because for steady wind I now also need to make a bts
                # file...
                caseName = makeName(
                    "", DOF, "U", URef, "RotatingWind", "", "PLExp", PLExp, "T", TMax
                )
            case "bts":
                case_dict["InflowFile|WindType"] = 3
                case_dict["InflowFile|FileName_BTS"] = f"Wind\\{windfile}"
            case "in":
                raise NotImplementedError(
                    "The OpenFAST toolbox does not have a function to read/write to .in files. Use a .dat file instead (you can just change the extension and everything should work ok)."
                )
            case "dat":
                f = FASTInputFile(templateDir / "Wind" / windfile)
                f["RandSeed1"] = randomSeed
                f["AnalysisTime"] = (
                    TMax + 300
                )  # Build in a bit of margin for the AnalysisTime.
                f["IECturbc"] = turbulenceIntencity
                f["URef"] = URef
                f["PLExp"] = PLExp
                newName = makeName(
                    "U",
                    URef,
                    "TI",
                    turbulenceIntencity,
                    "PLExp",
                    PLExp,
                    "seed",
                    randomSeed,
                    "T",
                    TMax,
                )
                newName += ".dat"
                f.write(templateDir / "Wind" / newName)

                case_dict["InflowFile|WindType"] = 3
                case_dict["InflowFile|FileName_BTS"] = (
                    f"Wind\\{newName.replace(".dat", ".bts")}"
                )
            case _:
                raise Exception(f"Couldn't match wind file extension {extension}.")

        case_dict["__name__"] = caseName

        params.append(case_dict)

    # Generate all the files in the work_dir
    case_gen.templateReplace(
        params,
        templateDir=str(templateDir),  # templateReplace uses strings rather than Paths.
        outputDir=str(outputDir),
        main_file=str(main_file),
        removeAllowed=True,  # Removes out, outb, ech, and sum files.
        removeRefSubFiles=True,  # I think this removes duplicate files.
        oneSimPerDir=True,
        dryRun=False,  # dryRun allows you to check the behaviour before finalizing the command.
    )

    # Do some final postprocessing that's not included in templateReplace.
    for param in params:
        caseName = param["__name__"]

        # If we use a bts wind file, it might be that we haven't run TurbSim on it yet.
        # That's because the .bts files are actually quite big so it is not good to
        # generate all possible .bts files and then copy them to all cases. So on a
        # per-case basis, we now generate a .bts file if needed.
        if param["InflowFile|WindType"] == 3:
            btsFile = Path(param["InflowFile|FileName_BTS"])
            if not btsFile.exists():
                datFile = outputDir / caseName / btsFile.with_suffix(".dat")
                cmd = f"src/TurbSim_x64.exe {str(datFile)}"
                subprocess.run(cmd, check=True)

        # templateReplace also copied the README file for each case, which is not
        # needed.
        readmeFile = Path(f"{outputDir}/{caseName}/README.md")
        readmeFile = outputDir / caseName / "README.md"
        readmeFile.unlink()


def makeName(*args) -> str:
    if not len(args) % 2 == 0:
        raise Exception("Must provide arguments in pairs of two.")

    nameParts = []
    for preFix, number in itertools.batched(args, 2):
        match number:
            case str() | int():
                numberName = str(number)
            case float():
                numberName = f"{number:.2g}".replace(".", "p")
            case _:
                raise TypeError("Unsupported type.")

        nameParts.append(preFix + numberName)
    return "_".join(nameParts)  # Remove the first "_".


if __name__ == "__main__":
    generate_cases(
        DOF="FullDOF",
        windfile="SimpleTurbulence.dat",
        randomSeeds=range(100, 110),
        turbulenceIntencities=[8],
        URefs=[15],
        PLExps=[0.07],
        TMaxs=[2100],
    )
