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