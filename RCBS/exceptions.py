class NotSingleAtomSelectionError(Exception):
    """
    Raised when the input selection is not a single atom
    """
    pass

class NotEqualListsLenghtError(Exception):
    """
        Raised when the 'sel_input' and the 'sel_types' list are not equal
    """
    def __init__(self):
        print("'sel_input' and the 'sel_types' list are not equal")
    pass

class NotEnoughAtomsSetectedError(Exception):
    """
        Raised when more than n atoms are required but not input.
    """
    def __init__(self):
        print("Number of selected atoms is not enough")
    pass

class  NotExistingMetalError(Exception):
    """
    Raised when the selected metal has no basis set (any of all the described) described.
    """

    def __init__(self):
        print('The metal has no basis set available')
    pass
