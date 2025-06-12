from os import path


def to_github(filename="", web=False, subdir=None):

    if web:
        directory = (
            "https://raw.githubusercontent.com/johnarban/wwt_interactives/main/images/"
        )
    else:
        directory = path.expanduser(f"~/github/wwt_interactives/images/")
    if subdir is not None:
        directory = path.join(directory, subdir)
    return path.join(directory, filename)

    #     return f'https://raw.githubusercontent.com/johnarban/wwt_interactives/main/images/{filename}'
    # return path.expanduser(f'~/github/wwt_interactives/images/{filename}')
