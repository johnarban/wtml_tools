DEBUG_LEVELS = {"DEBUG": 0, "INFO": 1}
DEBUG_LEVEL = 0


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
