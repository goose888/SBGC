"""

* File Name : auxiliary_lib.py

* Purpose : Lib stores approaches to execute general operations by calling 3rd-party sources / libs

* Creation Date : 18-07-2017

* Last Modified : Tue Jun 19 00:18:53 2018

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

def sweepdata(filename):
    """  sweep whole .csv file to check if all entries are in utf-8
    arguments: csv file name
    """
    data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID')
    for i in data.columns:
        row = 0
        for j in data[i]:
            row = row + 1    
            print "row = %d" % row
            try:
                str(j).decode('utf-8')
            except UnicodeError:
                print "string is not UTF-8, row = %d,  col = %s" % (row,i)
                print "value is ", j
                raw_input("Press Enter to continue...")

