from os import path


def to_github(filename="", web=False, subdir=None, repo="data_repo"):

    if web:
        directory = f"https://johnarban.github.io/data_repo/"
    else:
        directory = path.expanduser(f"~/github/data_repo/")
    if subdir is not None:
        directory = path.join(directory, subdir)
    return path.join(directory, filename)

    #     return f'https://raw.githubusercontent.com/johnarban/wwt_interactives/main/images/{filename}'
    # return path.expanduser(f'~/github/wwt_interactives/images/{filename}')
