
def get_html_table_header(columns, attributes):

    attributes_str = " ".join([ "%s=\"%s\"" %(x[0],x[1]) for x in attributes ])
    th_str = "<tr>%s</tr>" %"".join([ "<th>%s</th>" %x for x in columns ])
    return "<table %s>%s" %(attributes_str,th_str)
    #return "<table id=\"biana\" %s>%s" %(attributes_str,th_str)

def get_html_table_foot():
    return "</table>\n"

def append_html_table_values(values,rowIDs=None):

    if rowIDs is None:
        rowIDs = [ x for x in xrange(len(values)) ]

    data_str = "".join([ "<tr rowID=\"%s\">%s</tr>" %(rowIDs[current_row_index],"".join([ "<td>%s</td>" %str(current_column).replace("&","").replace("<","").replace(">","") for current_column in values[current_row_index] ])) for current_row_index in xrange(len(values)) ])
    return data_str

def get_html_table(columns,values,rowIDs=None,attributes=[],title=""):
    """
    "columns" must be a list with column headers
    """

    if rowIDs is None:
        rowIDs = [ x for x in xrange(len(values)) ]

    attributes_str = " ".join([ "%s=\"%s\"" %(x[0],x[1]) for x in attributes ])

    th_str = "<tr>%s</tr>" %"".join([ "<th>%s</th>" %x for x in columns ])
    data_str = "".join([ "<tr rowID=\"%s\">%s</tr>" %(rowIDs[current_row_index],"".join([ "<td>%s</td>" %str(current_column).replace("&","").replace("<","").replace(">","") for current_column in values[current_row_index] ])) for current_row_index in xrange(len(values)) ])
    return "<table %s>%s%s</table>\n" %(attributes_str,th_str,data_str)


def get_tabulated_table(values,columns=None):
    """
    """
    if columns is not None:
        to_print = ["\t".join(columns)]
    else:
        to_print = []

    to_print.extend(values)
    
    return "%s\n" %"\n".join( [ "\t".join(map(str,x)) for x in values ] )

