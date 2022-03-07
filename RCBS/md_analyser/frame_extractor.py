from ..snippets import check_folder


def frame_selector(
    frames_boolean,
    bins=10,
    frames_per_bin=1,
    bin_type="total_length",
    verbose=False,
    loop_tolerance=0.1,
):
    """
    DESCRIPTION:
        Function for randomly select frames from a boolean list (being the index the number of the frame).

    INPUT:
        - frames_boolean: list of boolean variables that indicate if the criteria is satisfied or not. The index of the list is the number of the frame.
        - bins:           number of partitions of the trajectory to properly distibute the selection
        - frames_per_bin: number of frames to select per bin
        - bin_type:       total_length | true_lenght --> total_length takes the whole set of frames to calculate the bins,
                            while true_length takes only the frames with True statements for calculating the bins
        - verbose:        True prints the frame number every time it is selected

    OUTPUT:
        - list of randomly selected frames
    """

    from random import choice

    selected = []

    frames = []
    for f in range(len(frames_boolean)):
        if frames_boolean[f] == True:
            frames.append(f + 1)
        elif frames_boolean[f] == False:
            pass

    if bin_type.lower() == "total_length":
        print("total_length")
        splits = int(round(len(frames_boolean) / bins, 0))

    elif bin_type.lower() == "true_length":
        print("true_length")
        splits = int(round(len(frames) / bins, 0))

    else:
        print(bin_type)

    for f in range(frames_per_bin):
        for b in range(bins):
            loopout = 0
            while True:
                if loopout + 1 == int(loop_tolerance * splits):
                    print(
                        "A frame between %s and %s has not been found. No frame will be saved for this bin."
                        % (splits * b, splits * (b + 1))
                    )
                    break

                else:
                    try:
                        frame = choice(range(splits * b, splits * (b + 1)))
                        if bin_type.lower() == "total_length":
                            if frame in frames:
                                selected.append(frame)
                                if verbose == True:
                                    print("The selected frame is:", frame)
                                break
                            else:
                                loopout += 1
                                continue

                        elif bin_type.lower() == "true_length":
                            selected.append(frames[frame])
                            break

                    except IndexError:
                        continue

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
