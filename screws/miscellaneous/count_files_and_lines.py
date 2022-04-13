



import os


def count_files_and_lines(start, files=0, lines=0, header=True, begin_start=None):
    """
    We use this function to count how many python files are in a path and how
    many lines of python code are in these files.

    The original code is provided by stack overflow user 'Bryce93' on page
    'https://stackoverflow.com/questions/38543709/count-lines-of-code-in-directory-using-python'.

    We have skipped
        1). files in folder '__PROGRAMS__'
        2). meaningless lines, i.e. of only spaces or '-'.
        3). files named with extension '.pyc'.

    Parameters
    ----------
    start : str
        The path we want to start tracing in.
    files :
    lines :
    header :
    begin_start :

    """
    if header:
        print("\t\t Line counting......\n")
        print('{:>5} |{:>8} |{:>10} | {:<20}'.format('NUM', 'ADDED', 'TOTAL', 'FILE'))
        print('{:->6}|{:->9}|{:->11}|{:->20}'.format('', '', '', ''))
    to_skip = [' ' * i + '\n' for i in range(30)]
    for thing in os.listdir(start):
        thing = os.path.join(start, thing)
        if os.path.isfile(thing):
            if thing.endswith('.py'):
                # we only look at python files.
                with open(thing, 'r') as f:
                    try:
                        FILES = f.readlines()[:]
                    except UnicodeDecodeError:
                        lines += 20
                    else:
                        newlines = len(FILES)
                        for entry in FILES:
                            if entry in to_skip:
                                newlines -= 1
                            else:
                                pass

                        lines += newlines
                        if begin_start is not None:
                            reldir_of_thing = '.' + thing.replace(begin_start, '')
                        else:
                            reldir_of_thing = '.' + thing.replace(start, '')

                        if len(reldir_of_thing) > 45:
                            reldir_of_thing = reldir_of_thing[:42] + '...'

                        files += 1
                        print('{:>5} |{:>8} |{:>10} | {:<20}'.format(int(files),
                                                                     newlines,
                                                                     lines,
                                                                     reldir_of_thing))

    for thing in os.listdir(start):
        # skip following ...
        if thing[-4:] == '.pyc':
            pass
        elif thing[:12] == '__contents__':
            pass
        elif thing[:4] == 'venv':
            pass
        else:
            thing = os.path.join(start, thing)
            if os.path.isdir(thing):
                files, lines = count_files_and_lines(thing,
                                                     files,
                                                     lines,
                                                     header=False,
                                                     begin_start=start)
    return files, lines


if __name__ == '__main__':
    import doctest
    doctest.testmod()

    count_files_and_lines('../../')
