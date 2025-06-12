DEBUG_LEVELS = {"DEBUG": 0, "STATS": 1, "INFO": 2, "ERROR": 3, "ALWAYS": 4}
DEBUG_LEVEL = 3


def set_debug_level(level):
    global DEBUG_LEVEL
    if isinstance(level, str):
        DEBUG_LEVEL = DEBUG_LEVELS[level]
    else:
        DEBUG_LEVEL = level


def log(arg: str, level=0):
    # check if level is a string
    if isinstance(level, str):
        level = DEBUG_LEVELS[level]

    if level >= DEBUG_LEVEL:
        print(arg)
