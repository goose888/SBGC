"""

* File Name : auxiliary_lib.py

* Purpose : Lib stores approaches to execute general operations by calling 3rd-party sources / libs

* Creation Date : 18-07-2017

* Last Modified : Wed 19 Jul 2017 03:01:54 AM EDT

* Created By : Shijie Shu

"""

import subprocess as sp

def exec_shell(command, path_of_script, args):
    """ Execute command through linux shell
    Input:
        command --- The command being used under linux shell
        path_of_script --- The name and path of the script being called by the command
        args --- The list of arguments being used in the script
    Output:
        status --- Return the standard output of the command
    """

    # Build subprocess command
    cmd = [command, path_of_script] + args
    # check_output will run the command and store to result
    status = sp.check_output(cmd, universal_newlines=True)

    return status

def rep_sp_by_sym(data, sym):
    """ Replace the space in the column name by the specified symbol
    Input:
        data --- The pandas data frame with column names already spcified
    Output:
        cols --- Return the new column name after replacement
    """

    cols = data.columns
    cols = cols.map(lambda x: x.replace(' ', sym) if isinstance(x, (str, unicode)) else x)

    return cols


