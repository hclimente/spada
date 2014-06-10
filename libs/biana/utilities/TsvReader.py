
from FormattedFileProcessor import FormattedFileProcessor

class TsvReader(FormattedFileProcessor):
    """
	Read/process TSV (tab seperated) formatted files
    """
    def __init__(self, input_file_name):
	FormattedFileProcessor.__init__(self, input_file_name=input_file_name, input_type="tsv")
	return

    def process(self, out_method, fields_to_include, overwrite_keys, keys_to_include):
	"""
	    Read and process an input file line by line. If out_method is None a dictionary storing read lines are returned.
	    out_method: method to output columns in current line on the fly in tsv format
	    fields_to_include: columns that would be included in the dictionary or processed with the function
	    overwrite_keys: allows overwriting keys (in case of entries with duplicate primary column values in the file) returning values as list in the dictionary
			    if False, returns list of values as list in the dictionary (each list element corresponding to value list of distinct entries)
	    keys_to_include: use only lines whose value of the first column is inside this set. Set None for including all the lines in the file
	"""
	file = open(self.input_file_name)
	line = file.readline()
	cols = [ c.lower() for c in line.strip().split('\t') ]
	if fields_to_include is None:
	    first_column = cols[0]
	else:
	    fields_to_include = [ f.lower() for f in fields_to_include ]
	    first_column = fields_to_include[0]
	columns = dict(zip(cols, range(len(cols))))
	#print "Start:", columns
	id_to_value = {}
	i=0
	line_prev = line
	line = file.readline()
	vals = line.strip().split('\t')
	#print vals
	while line:
	    try:
		vals = line.strip().split('\t')
		id = vals[columns[first_column]]
		if keys_to_include is None or id in keys_to_include:
		    if out_method is None:
			if fields_to_include is None:
			    if overwrite_keys:
				id_to_value[id] = vals
			    else:
				id_to_value.setdefault(id, []).append(vals)
			else:
			    if overwrite_keys:
				id_to_value[id] = [ vals[columns[f]] for f in fields_to_include]
			    else:
				id_to_value.setdefault(id, []).append( [vals[columns[f]] for f in fields_to_include] )
		    else:
			out_method("%s\n" % "\t".join([ vals[columns[f]] for f in fields_to_include ]))
	    except: #Exception, e: 
	    	#print "In: ", __file__, e, vals
		print line_prev, line 
		import traceback
		traceback.print_exc()
	    i+=1
	    #if i > 20:
	    #	break
	    line_prev = line
	    line = file.readline()
	file.close()
	#print columns, "\n", id_to_value.items()[0]
	if out_method is None:
	    if fields_to_include is not None:
		cols2 = []
		for c in cols:
		    if c in fields_to_include:
			cols2.append(c)
		columns = dict(zip(cols2, range(len(cols2))))
	    #print "End:", columns
	    return columns, id_to_value
	else:
	    return


