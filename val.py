import os

def is_valid_file(filename):
    if(os.path.exists(filename)):
        return True
    return False
