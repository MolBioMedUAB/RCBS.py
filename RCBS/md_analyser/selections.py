from ..exceptions import NotEqualListsLenghtError


def selection(u, sel_input, sel_type=None):
    """
    DESCRIPTION
        This function takes an input number or name and type of selection (atom or residue or none) and transforms it into an MDAnalysis selection. The input can also be a list of numbers or name and it can combine both types of input by using a list of selection types (of the same lenght).

    INPUT
        - u: MDAnalysis' universe
        - sel_input: int, str, list or tuple of int, or list or tuple of str
        - sel_type: type of selection input as sel_input
            at_num   -> atom number
            at_name  -> atom name
            res_num  -> residue number
            res_name -> residue name
            none     -> pipes the selection command directly

    OUTPUT
        - MDAnalysis selection as AtomGroup
    """

    if sel_type == None:
        sel_string = sel_input

    elif isinstance(sel_input, list):
        if isinstance(sel_type, str):
            if sel_type.lower() == "at_num":
                sel_string = "bynum " + " or bynum ".join(str(at) for at in sel_input)
            elif sel_type.lower() == "at_name":
                sel_string = "name " + " or name ".join(sel_input)
            elif sel_type.lower() == "res_num":
                sel_string = "resid " + " or resid ".join(sel_input)
            elif sel_type.lower() == "res_name":
                sel_string = "resname " + " or resname ".join(sel_input)
            else:
                sel_string = " or ".join(sel_input)

        elif isinstance(sel_type, list):
            if len(sel_input) != len(sel_type):
                raise NotEqualListsLenghtError

            else:
                sel_string = ""
                for inp, typ in zip(sel_input, sel_type):
                    if sel_string == "":
                        pass
                    else:
                        sel_string = "".join([sel_string, " or"])

                    if typ.lower() == "at_num":
                        sel_string = " ".join([sel_string, "bynum", str(inp)])
                    elif typ.lower() == "at_name":
                        sel_string = " ".join([sel_string, "name", str(inp)])
                    elif typ.lower() == "res_num":
                        sel_string = " ".join([sel_string, "resid", str(inp)])
                    elif typ.lower() == "res_name":
                        sel_string = " ".join([sel_string, "resname", str(inp)])
                    else:
                        sel_string = " ".join([sel_string, str(inp)])

    elif isinstance(sel_input, (str, int)):
        if sel_type.lower() == "at_num":
            sel_string = " ".join(["bynum", str(sel_input)])
        elif sel_type.lower() == "at_name":
            sel_string = " ".join(["name", str(sel_input)])
        elif sel_type.lower() == "res_num":
            sel_string = " ".join(["resid", str(sel_input)])
        elif sel_type.lower() == "res_name":
            sel_string = " ".join(["resname", str(sel_input)])

    return u.select_atoms(sel_string)
