


from screws.decorators.accepts import accepts
import re



@accepts(str)
def check_no_splcharacter(test):
    """Check the string ``test`` have no special character listed below.

    :param test:
    :return:
    """
    # noinspection RegExpRedundantEscape
    string_check = re.compile('[@!#$%^&*()<>?/\|}{~:]')
    if string_check.search(test) is None:
        return True
    else:
        return False