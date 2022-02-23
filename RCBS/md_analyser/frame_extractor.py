from ..snippets import check_folder


def frame_selector(frames_boolean, bins=10, frames_per_bin=1):
    """
    DESCRIPTION:
        Function for randomly select frames from a boolean list (being the index the number of the frame).

    INPUT:
        - frames_boolean: list of boolean variables that indicate if the criteria is satisfied or not. The index of the list is the number of the frame.
        - bins:           number of partitions of the trajectory to properly distibute the selection
        - frames_per_bin: number of frames to select per bin

    OUTPUT:
        - list of randomly selected frames
    """

    from random import choice

    selected = []
    splits = int(round(len(frames_boolean) / bins, 0))

    for f in range(frames_per_bin):
        for b in range(bins):
            while True:
                frame = choice(range(splits * b, splits * (b + 1)))
                try:
                    if frames_boolean[frame] == True:
                        break
                    elif frames_boolean[frame] == False:
                        continue

                except IndexError:
                    continue

            selected.append(frame + 1)

    selected = sorted(selected)

    return selected


def trajectory_frame_extractor(
    u, frame, name="", folder="", outformat=".pdb", solvent=True
):
    """
    DESCRIPTION:
        Function for saving a frame from a trajectory loaded as a MDAnalysis Universe

    INPUTS:
        u       --> MDAnalysis Universe
                --> It has to be a trajectory.
        frame   --> it can be a single frame (input as number or string of number) or a list of frames.
        name    --> name for the output file. If more than one frame is input, outputs will be secuencially named
        folder  --> folder name where pdbs have to be saved.
                --> if the folder does not exist, it will be created.
                --> by default, no folder is used
        format  --> format for the output file
                --> by default, output will be saved as pdb. All formats accepted by the write module of MDAnalysis are accepted.

    """

    structure_saver = lambda sel, file_name, outformat: sel.write(file_name + outformat)

    if folder != "":
        check_folder(folder)

    elif folder == "":
        folder = "."

    if isinstance(frame, list):
        f_ = 0
        for f in frame:
            f_ += 1
            u.trajectory[int(f)]
            if solvent == True:
                sel = u.select_atoms("all")
            elif solvent == False:
                sel = u.select_atoms(
                    "all and not (resname WAT or resname NA+ or resname Na+ or resname Cl- or resname CL-)"
                )

            if name == "":
                file_name = folder + "/" + str(f)
            else:
                file_name = folder + "/" + name + "_" + str(f_)

            structure_saver(sel, file_name, outformat)

    if isinstance(frame, (str, int)):
        u.trajectory[int(frame)]
        if solvent == True:
            sel = u.select_atoms("all")
        elif solvent == False:
            sel = u.select_atoms(
                "all and not (resname WAT or resname NA+ or resname Na+ or resname Cl- or resname CL-)"
            )

        if name == "":
            file_name = folder + str(frame)
        else:
            file_name = folder + name

        structure_saver(sel, file_name, outformat)
