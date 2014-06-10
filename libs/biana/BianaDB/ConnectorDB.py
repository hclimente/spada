"""
    BIANA: Biologic Interactions and Network Analysis
    Copyright (C) 2009  Javier Garcia-Garcia, Emre Guney, Baldo Oliva

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

"""
 File       : ConnectorDB.py, derived from file PianaDB.py from the PIANA project
 Author     : Javier Garcia & Emre Guney, based on R. Aragues & D. Jaeggi script
 Creation   : 2003
 Contents   : class for establishing conexions to BianaDatabase and handling inserts and selects
 Called from: BianaDBaccess.py
=======================================================================================================

This a generalization of mysql commands

To see how this class is used, look at any method in BianaDBaccess.py

"""

#import biana.ext.MySQLdb as MySQLdb
import mysql.connector as db_connector
import sys

DEBUG_BUFFER_INSERT_SINGLE = False # Set True to control queries, will insert each query seperately
DEBUG_PRINT_INSERT_QUERY = False  # Set True to control queries, will print insert_db_content queries
DEBUG_CHECKING_TABLES = False

class DB(object):
    """
    Class for establishing conexions to pianaDB and handling inserts and selects
    """

    def __init__(self, dbname=None, dbhost=None, dbuser=None, dbpassword=None, dbport=None, dbsocket=None, buffer=True, lock_tables=False):
        """
        "dbname" is the database name to which you want to connect to (required)
 
        "dbhost" is the machine with the mysql server that holds the piana database (required)

        "dbuser" is the mysql user (not required in most systems)

        "dbpassword" is the mysql password (not required in most systems)

        "dbport" is the mysql port (not required in most systems)

        "buffer" determines if an insert buffer is going to be used. Default is True. Alert! DB Object must be closed before in order to empty the buffer!
        Buffer is used due to performance when parsing data

        "lock_tables" is used to lock the used tables
        """

        self.dbname = dbname
        self.dbhost = dbhost
        self.dbuser = dbuser
        self.dbpassword = dbpassword
        self.dbport = dbport
        self.dbsocket = dbsocket

        # init database connection (different connection parameters depending on user preferences...)
        if dbsocket is not None or dbsocket=="":
            if not dbport:
                if not dbuser is None and not dbpassword is None:
                    self.db = db_connector.connect(host= dbhost, user=dbuser, passwd=dbpassword, unix_socket=self.dbsocket)
                elif not dbuser is None:
                    self.db = db_connector.connect(host= dbhost, user=dbuser, unix_socket=self.dbsocket)
                else:
                    self.db = db_connector.connect(host= dbhost, unix_socket=self.dbsocket)
            else:
                if not dbuser is None and not dbpassword is None:
                    self.db = db_connector.connect(host= dbhost, user=dbuser, passwd=dbpassword, port=dbport, unix_socket=self.dbsocket)
                elif not dbuser is None:
                    self.db = db_connector.connect(host= dbhost, user=dbuser, port=dbport, unix_socket=self.dbsocket)
                else:
                    self.db = db_connector.connect(host= dbhost, port=dbport, unix_socket = self.dbsocket)
        else:
            if not dbport:
                if not dbuser is None and not dbpassword is None:
                    self.db = db_connector.connect(host= dbhost, user=dbuser, passwd=dbpassword)
                elif not dbuser is None:
                    self.db = db_connector.connect(host= dbhost, user=dbuser)
                else:
                    self.db = db_connector.connect(host= dbhost)
            else:
                if not dbuser is None and not dbpassword is None:
                    self.db = db_connector.connect(host= dbhost, user=dbuser, passwd=dbpassword, port=dbport)
                elif not dbuser is None:
                    self.db = db_connector.connect(host= dbhost, user=dbuser, port=dbport)
                else:
                    self.db = db_connector.connect(host= dbhost, port=dbport)


        # Not necessary to do autocommit as the Engine selected for tables is MyISAM (MyISAM does not accept commit)
        # If anytime this changes, it will be necessary to activate autocommit or to do commit at each parser
        #self.db.autocommit(1)

        self.cursor = self.db.cursor()

        self.dbmaxpacket = self._get_max_packet()
        self.lock_frequency = 100 #20000
        self.current_lock_num = 0
        self.is_locked = False
        self.locked_tables = set()
        self.locked_tables_win = set() # The same as self.locked_tables, but corrected for windows (as table names are case insensitive and they usually are lowercased)
        self.lock_tables = lock_tables

        # Set dbmaxpacket
        #self.insert_db_content("SET SESSION max_allowed_packet=16777216")
        #self.insert_db_content("SET SESSION max_allowed_packet=4777216")
        #self.insert_db_content("SET SESSION max_allowed_packet=1000000")
        #self.insert_db_content("SET SESSION max_allowed_packet=1047552")

        if buffer is True and self.dbmaxpacket is not None:
            self.uses_buffer = True
            self.insert_buffer = Buffer(self.dbmaxpacket,self)
            self.autoincrement_values = {}
        else:
            self.uses_buffer = False
            self.insert_buffer = None

        if( dbname is not None ):
            #self.select_db_content("use "+dbname)
	    self.use_database(dbname)

        self.table_names = None

        return


    def use_database(self, database_name):
        """
        Specifies which database to use
        """
        self.db.database = database_name
	self.dbname = database_name


    # --
    # methods needed for using pickle with piana objects
    # -- 

    def __getstate__(self):

        odict = self.__dict__.copy() # copy the dict since we are going to change it
        del odict['db']              # remove conexion to MySQL: attribute self.db cannot be pickled
        return odict

    def __setstate__(self, dict):

        self.__dict__ = dict          # recover previous dictionary 
        dict.__class__.__init__(dict) # recover conexion to MySQL by calling init
        
    def __getnewargs__(self):

        return (self.dbname, self.dbhost, self.dbuser, self.dbpassword)  # returns arguments that init will get
                                                                         # when it is called after the pickle.load

    def __str__(self):

        return "Connection to database %s on %s as %s" %(self.dbname, self.dbhost, self.dbuser)


    def close(self):
        """
        Closes the connection with the database
        """
        self._unlock_tables()
        self.cursor.close()
        self.db.close()

    def check_consistency_with_given_source_version(self, source_code_version):
        #print self._get_source_code_version_in_db(), source_code_version
        if self._get_source_code_version_in_db() is None: # newly created database
            self._set_source_code_version_in_db(source_code_version)
            return True
        elif self._get_source_code_version_in_db().lower() == source_code_version.lower(): # check code version
            return True
        return False

    # --
    # handling inserts and selects through a generalized class
    # --

    def insert_db_content(self, sql_query, answer_mode = None):
        """
        Inserts values into a piana database (connection was established in self.db)

        Depending on argument "answer_mode", different things are returned.

        This method is called by PianaDBaccess to then process the results and return them to the user

        "slq_query" is normally obtained through classes implemented in PianaInsertSQL.py, which have a method get_sqlquery
        that creates the sql query needed to retrieve the searched value.

        "answer_mode" can be one of the following:

        - None: nothing is returned
        - 'last_id' : last id inserted is returned
        - 'num_updated' : number of rows that were updated (used by UPDATE statements...)

          --> 'last_id' mode only works for those tables that have an auto_increment ID!!!!!!! Will not work with primary keys that are not
              auto_increment. Currently, following tables have auto_increment ids: protein (ID=proteinPiana) and interaction (ID=interactionPiana)

        "buffer" indicates if this sql_query has to be inserted or not. It can be:
        - None: sql_query will be treated as a sql_query and it will be executed (if possible)
        - dictionary: the insert has been inserted into the dictionary buffer. So, the sql_query will not be executed. (So, sql_query can be None too. It will be ignored.

        """

        if sql_query is None:
            return None

        # Checks lock frequency
        #self._check_lock_frequency()  # ENTERS INTO A LOOP
   
        if isinstance(sql_query,list):
            for actual_query in sql_query:
                if DEBUG_PRINT_INSERT_QUERY:
                    sys.stderr.write(actual_query+"\n")
                try:
                    if( actual_query != "" ):
                        self.cursor.execute(actual_query)
                except Exception, inst:
                    sys.stderr.write("Attention: this query was not executed due to a mysql exception: <<%s>>\n" %(actual_query))
                    sys.stderr.write("           Error Reported: %s\n" %(inst))
                    answer_mode = None
                    raise ValueError(inst)
            return None

                
        if sql_query is not None:
            if DEBUG_PRINT_INSERT_QUERY:
                sys.stderr.write("%s\n" %sql_query)
            try:
                #sys.stderr.write("Sending query to MySQL server (%s)...\n" %(sql_query) )
                # In order to avoid any insertion
                #pass
                self.cursor.execute(sql_query)
            except Exception, inst:
                sys.stderr.write("Attention: this query was not executed due to a mysql exception: <<%s>>\n" %(sql_query))
                sys.stderr.write("           Error Reported: %s\n" %(inst))
                raise ValueError(inst)
                answer_mode = None
        #else:
        #    raise ValueError("Trying to execute an empty insertion sqlquery")

        if answer_mode == "last_id":
            # in case mode "last_id" was chosen, do a select to the database to obtain the last autoincrement id inserted
            
            self.cursor.execute("""Select LAST_INSERT_ID()""")
            answer = self.cursor.fetchall()

            if answer:
                return answer[0][0]
            else:
                return None
        
        elif answer_mode == "num_updated":
            # in case mode "num_updated" was chosen, the number of columns that were updated will be in the answer of the sql query
            answer = self.cursor.fetchall()
            if not answer:
                return None
            else:
                return answer[0]
        
        else:
            return None
        # END OF else: (elif answer_mode == "num_updated":)

        sys.stderr.write("Query executed!\n")

        
    def select_db_content(self, sql_query= None, answer_mode="single", remove_duplicates="yes", number_of_selected_elems= 1):
        """
:        Returns content from a piana database (connection was established in self.db)
        
        "slq_query" is normally obtained through classes implemented in PianaSelectSQL.py, which have a method get_sqlquery
        that creates the sql query needed to retrieve the searched value
        
        "answer_mode" is used to let the user choose what kind of answer he is expecting from the sql query
        answer_mode can take values (default is "single"):
        - "single": just one element (if nothing is found, returns None)
        - "list": a list of single elements (if nothing is found, returns empty list [])
        - "raw":  the raw result from the sql query (a matrix)
        
        "remove_duplicates" is used to let the user decide if the returning list can contain duplicates or not
        (only used when answer_mode="list")
        (this is useful to remove similar entries from one query, for example same uniprot accession numbers returned
        that are actually the same under different swissAccession source DB )
        
        - "yes" will return a list where all elements are unique
        - "no" will return the complete list returned by the sql query
        
        "number_of_selected_elems" sets the number of elements being selected in the sql query.
        - If 1, only the string of the element is returned. 
        - If >1, each group of elements is represented by a tuple (elem1, elem2, ...)
        - if you want to use a number_of_selected_elems higher than 3, you have to modify the code below
        to create keys of the number of elements you wish to select
        """

        if sql_query is not None:
            try:
                self.cursor.execute(sql_query)
                answer = self.cursor.fetchall()
            except Exception, inst:
                sys.stderr.write("Attention: this query was not executed due to a mysql exception: <<%s>>\n" %(sql_query))
                sys.stderr.write("           Error Reported: %s\n" %(inst))
                answer = None
                #print self._get_max_packet()
                raise ValueError(inst)
        else:
            raise ValueError("Trying to execute an empty selection sqlquery")
            
        # when fetchall doesn't find anything, it returns an empty object.
        # if something was found, transform to answer_mode requested by user
        
        # TO DO!!!! Speed up things by creating a dictionary with those values that have been already inserted...
        #           that will be much faster than checking on a list!!!
        
        if answer:
            #sys.stderr.write("Answer: %s Answer[0]: %s Answer[0][0]: %s\n" %(answer, answer[0], answer[0][0]) )
            # something was found: each answer_mode returns a different "something"
            if answer_mode == "single":
                if number_of_selected_elems == 1:
                    return answer[0][0]
                else:
                    return answer[0]
            
            elif answer_mode == "list":
                # Optimization? I don't know... Done by Javi. It is not necessary to sort to form the key (MySQL always will return the result in the same order...)
                # Furthermore, in many tables data is sorted in some ways (for example, in interactionTable proteinPianaA and proteinPianaB are sorted...
                # Take into account this!

                if number_of_selected_elems == 1:
                    temporal_list = [element[0] for element in answer]
                    if remove_duplicates=="yes":
                        return list(set(temporal_list))
                    else:
                        return temporal_list
                    
                else:
                    # if there are more than one selected elements:
                    #   -> 
                    #      
                    temporal_list = list(answer)

                    if remove_duplicates=="yes":
                        #return dict(map(None,temporal_list,[None])).keys()
                        return dict(map(None,temporal_list,[])).keys()
                    else:
                        return temporal_list

            elif answer_mode == "raw":
                return answer
            else:
                raise ValueError("answer_mode value is not correct")
        # END OF if answer:
        
        else:
            # nothing was found: each answer_mode returns a different "nothing"
            if answer_mode == "single":
                return None
            elif answer_mode == "list":
                return []
            elif answer_mode == "raw":
                return answer
            else:
                raise ValueError("answer_mode value is not correct")
        # END OF else: (if answer:)


    def add_autoincrement_columns(self, table, attribute):
        
        if self.uses_buffer is False:
            raise ValueError("Cannot add an autoincrement column if buffer is not being used")
        
        if self.autoincrement_values.has_key((table,attribute)):
            raise ValueError("Trying to add twice the same autoincrement column")

        self.autoincrement_values[(table,attribute)] = None


    def get_next_autoincrement(self, table, attribute ):

        self._check_locked_table(table)

        if self.is_locked is False:
            raise ValueError("Trying to get an autoincrement value without locking the table")

        if self.autoincrement_values[(table,attribute)] is None:
            self.autoincrement_values[(table,attribute)] = self._get_current_autoincrement(table, attribute)
        
        self.autoincrement_values[(table,attribute)] += 1
        
        return self.autoincrement_values[(table,attribute)]

    def _get_last_stable_autoincrement(self, table, attribute):
        return self.select_db_content(self._get_select_sql_query(tables=["BianaDatabase"], columns=["last_"+attribute]))

    def _get_current_autoincrement(self, table, attribute ):

        if self.is_locked is False:
            raise ValueError("Trying to get an autoincrement value without locking the table")

        if self.autoincrement_values[(table,attribute)] is None:
            #max_value = self._get_max_value(table,attribute)
            #if max_value is None:
            #    return 0
            #else:
            #    return max_value
            #max_value = self.select_db_content("SELECT %s FROM %s" % (attribute, table), answer_mode="single" )
            max_value = self._get_last_stable_autoincrement( table = table,
                                                             attribute = attribute )
            #print max_value
            return max_value
        else:
            return self.autoincrement_values[(table,attribute)]


    def set_lock_tables(self, value):
        self.lock_tables = value
        
    def _check_locked_table(self, table):
    	if self.lock_tables is True:
            if table.lower() not in self.locked_tables_win:
                self.locked_tables.add(table)
                self.locked_tables_win.add(table.lower())
                self._unlock_tables()
                self._lock_tables()

    def _get_source_code_version_in_db(self):
        """
        Fetches value in source_code_version field in BianaDatabase table
        """
        table = "BianaDatabase"
        column = "source_code_version"
    	self._check_locked_table(table)
        return self.select_db_content( "SELECT %s FROM %s" % (column, table), answer_mode="single" )

    def _set_source_code_version_in_db(self, source_code_version):
        """
        Sets value in source_code_version field in BianaDatabase table
        """
        table = "BianaDatabase"
        column = "source_code_version"
    	self._check_locked_table(table)
        return self.insert_db_content( "UPDATE %s set %s=\"%s\"" % (table, column, source_code_version) )

    def _get_max_value(self, table, column):
        """
        OBSOLETE - last_ fields in BianaDatabase is used now
        """

    	self._check_locked_table(table)
        return self.select_db_content( "SELECT MAX(%s) FROM %s" %(column,table), answer_mode="single" )

    def _get_max_packet(self):
        """
        Returns the maximum packet size that the sql server accepts
        """

        query = "show variables like \"max_allowed_packet\""

        maxpacket_res = self.select_db_content(sql_query= query, answer_mode="raw", remove_duplicates="yes", number_of_selected_elems= 1)

        if maxpacket_res:
            return maxpacket_res[0][1]
        else:
            sys.stderr.write("Unknown database maximum packet size")
            return None


    def check_database(self, database, ignore_primary_keys=False, verbose=False):
        """
        Method for checking database

        It checks if all tables are existing

        It does not drop any table

        For the moment, it does not check default value of table fields
        """
        
        #return

        existing_tables_str = set(self._get_table_names())
        
        if( sys.platform.lower()=="win32" or sys.platform.lower()=="win64" ):
            existing_tables_str = set([ x.lower() for x in existing_tables_str ])

        checked_tables_str = set()

        database_tables = database.get_tables()

        for current_table_obj in database_tables:
            current_table_str = current_table_obj.get_table_name()

            if( sys.platform.lower()=="win32" or sys.platform.lower()=="win64" ):
                current_table_str = current_table_str.lower()

            checked_tables_str.add(current_table_str)
            if current_table_str not in existing_tables_str:
                # Table does not exist. Create it
                self.insert_db_content( sql_query = current_table_obj.create_mysql_query(ignore_primary_key=ignore_primary_keys) )
            else:
                # Table exists. Check column names & types
                existing_columns = self.select_db_content( sql_query = "DESC %s" % current_table_str,
                                                           answer_mode = "raw" )
                # It should be:
                # 0.- field name
                # 1.- Type
                # 2.- Null
                # 3.- Key
                ## differantiating between table column (info comes from database) and table field (info comes from database.py)
                column_dictionary = {}
                for current_column in existing_columns:
                    column_name = current_column[0].lower()
                    ## !!mysql (at least 5.0.22) automatically converts int->int(10), smallint->smallint(5), float(4)->float and text(xxx)->text, bool->tinyint(1)  
                    column_type = current_column[1].lower().replace("smallint(5)", "smallint").replace("int(10)", "int").replace("float(4)", "float").replace("integer", "int").replace("'","\"").strip() 
                    if current_column[2].startswith("NO"): #find("NO") != -1:
                        column_null = False
                    elif current_column[2].startswith("YES"):
                        column_null = True
                    else:
                        print "Warning: NULL column not recognized: ", current_column[2]
                    column_extra = current_column[5].lower()
                    if column_extra != "":
                        column_type += " " + column_extra
                    column_dictionary[column_name] = (column_type, column_null)
                table_fields = current_table_obj.get_fields()
                for current_field in table_fields:
                    field_name = current_field.field_name.lower()
                    field_type = current_field.data_type.lower().replace("integer", "int").replace("float(4)", "float").replace("bool","tinyint(1)").replace("'","\"").replace("default 0","").strip()
                    field_null = current_field.null
                    ## check whether field name is existing
                    if column_dictionary.has_key(field_name):
                        ## check whether field data is consistent
                        (column_type, column_null) = column_dictionary[field_name]
                        if field_type != column_type or field_null != column_null:
                            if field_null == column_null and column_type == "text":
                                continue
                            #if verbose:
                            #    print "NEW: ",field_type, column_type, field_null, column_null
                            if field_null:
                                null_str = "NULL"
                            else:
                                null_str = "NOT NULL"
                            self._check_locked_table( current_table_str )    
                            query_str = "ALTER TABLE %s MODIFY %s %s %s" % (current_table_str, current_field.field_name, field_type, null_str) 
                            self.insert_db_content( sql_query = query_str )
                            #print field_type, column_type, field_null, column_null, query_str
                            if verbose:
                                sys.stderr.write("%s\n" %query_str)
                    else:
                        if field_null:
                            null_str = "NULL"
                        else:
                            null_str = "NOT NULL"
                        self._check_locked_table( current_table_str )
                        query_str = "ALTER TABLE %s ADD %s %s %s" % (current_table_str, current_field.field_name, field_type, null_str) 
                        self.insert_db_content( sql_query = query_str )
                        #print field_type, column_type, field_null, column_null, query_str
                        if verbose:
                            #print "NEW: ",field_type, column_type, field_null, column_null
                            sys.stderr.write("%s\n" %query_str)


        if DEBUG_CHECKING_TABLES:
            checked_tables_str = existing_tables_str - checked_tables_str
            if len(checked_tables_str)>0:
                for db_str in checked_tables_str:
                    if not db_str.lower().startswith("userentityunification_protocol_") and not db_str.lower().startswith("key_attribute_"):
                        print "Warning: Following tables are not encoded in source: ",db_str
        return


    def GetTableNames(self):
        """
        Generates the SQL statement that gets all database table names
        """
        
        return """ SHOW TABLES """

    
    def _get_select_sql_query(self,tables,columns=None,fixed_conditions=None,join_conditions=None, group_conditions=None, distinct_columns=False):
        """
        Generates a general select sql statement

        "tables" is a list or tuple of tables where the value/s must be searched. If the elements of the list or tuple are tuples of length 2, the format taken will be the following:

                      (table_name or table_object, alias to the table)

        "columns" is a list or tuple of the columns searched (columns must be preceeded by the table where they are searched). If it is None, all values will be selected

        "fixed_conditions" is a list or tuple of tuples with the following format: (column,type,restriction_value)

        "join_conditions" is a list or tuple of tuples with the following format: (column,type,column) to restrict the selection to the joint

        "type" can be "=",">","<",...

        "group_conditions" is a list of columns where it must be grouped 

        It returns the sql query.
        """
       

        if( fixed_conditions is None or not len(fixed_conditions) ):
            fixed_conditions_sql = ""
        else:
            #fixed_conditions_sql = " AND ".join(["%s %s \"%s\"" %(x[0],x[1],x[2]) for x in fixed_conditions])
            fixed_conditions_sql = " AND ".join( [self._get_cond_str(x) for x in fixed_conditions] )

        if( join_conditions is None or not len(join_conditions) ):
            join_conditions_sql = ""
        else:
            join_conditions_sql = " AND ".join(["%s %s %s" %(x[0],x[1],x[2]) for x in join_conditions])
            if fixed_conditions_sql != "":
                join_conditions_sql = " AND %s" %(join_conditions_sql)

        if( join_conditions or fixed_conditions ):
            where_sql = " WHERE "
        else:
            where_sql = ""


        if columns is None:
            columns_sql = "*"
        else:
            columns_list = []
            for current_column in columns:
                if( isinstance(current_column, tuple) ):
                    columns_list.append("%s AS %s" %(current_column))
                else:
                    columns_list.append(current_column)
            columns_sql = ",".join(columns_list)

        # tranform table objects to table name strings
        #tables = [ "%s" %(x) for x in tables ]
        # Check if the tables are locked in case it is necessary
        tables_list = []
        for actual_table in tables:
            if( isinstance(actual_table,tuple) ):
                tables_list.append("%s AS %s " %(actual_table[0],actual_table[1]) )
                self._check_locked_table(str(actual_table[1]))
            else:
                tables_list.append("%s" %actual_table)
                self._check_locked_table(str(actual_table))
                		
        if group_conditions is not None and len(group_conditions)>0:
            group_conditions_sql = "GROUP BY %s" %(",".join(group_conditions))
        else:
            group_conditions_sql = ""

        if distinct_columns:
            distinct_str = "DISTINCT"
        else:
            distinct_str = ""
            
        return """SELECT %s %s FROM %s %s %s %s %s""" %(distinct_str,
                                                        columns_sql,
                                                        ",".join(tables_list),
                                                        where_sql,
                                                        fixed_conditions_sql,
                                                        join_conditions_sql,
                                                        group_conditions_sql)
    
    def _get_delete_sql_query(self, table, fixed_conditions=None):
        """
        Generates a general delete sql statement
        
        "table" is a table name or table object

        "fixed_conditions" is a list or tuple of tuples with the following format: (column,type,restriction_value)

        "type" can be "=",">","<",...

        It returns the sql query.
        """
       
        self._check_locked_table(table)

        if( fixed_conditions is None or not len(fixed_conditions) ):
            fixed_conditions_sql = ""
        else:
            fixed_conditions_sql = " AND ".join( [self._get_cond_str(x) for x in fixed_conditions] )

        if fixed_conditions:
            where_sql = " WHERE "
        else:
            where_sql = ""

        return """DELETE FROM %s %s %s""" %(table,
                                            where_sql,
                                            fixed_conditions_sql)
                                                  

    def _get_union_queries(self, queries):
        """
        "queries" is the list of queries to make the union
        """

        return " UNION ".join(queries)

    
    def _get_cond_str(self, conditions):
        
        if len(conditions) == 4:
            return "%s %s %s" %(conditions[0],conditions[1],conditions[2])
        else:
            return "%s %s \"%s\"" %(conditions[0],conditions[1],conditions[2])


    ####
    # INSERT RELATED METHODS
    ####

    def _get_nested_insert_sql_query(self, table, columns, subquery):

        self._check_locked_table(table)
        
        return """INSERT INTO %s (%s) (%s)""" %(table, ",".join(columns), subquery)


    def _get_insert_sql_query(self,table,column_values,special_column_values=[],use_buffer=True, max_elements_in_buffer=None, on_duplicate_key="IGNORE"):
        """
        Generates a general sql statement
        
        "column_values" must be a tuple of tuples with the following format: (column, value)
    
        It returns the sql query. It does not use the buffer!
        
        By default, if the buffer is active, it will use the buffer, unless the "use_buffer" atributte is set to False

        It can return a list of queries instead a single query in the case the buffer is being used

        "max_elements_in_buffer" is only used when buffer is used. It is not mandatory

        "on_duplicate_key" specifies the query to be used when duplicate keys exists
        
        """

        # tranform table objects to table name
        table = "%s" %(table)

        self._check_locked_table(table)
        
        columns = []
        values = []
        for x in column_values:
            columns.append(x[0])
            if x[1] is None:
                raise ValueError("Trying to insert a None in table %s" %(table))
            if isinstance(x[1],unicode):
                # Transform to ascii
                value = x[1].encode('ascii','replace')
            else:
                value = str(x[1])
            values.append("\"%s\"" %value.replace('\\','\\\\').replace('"','\\"'))

        #columns, values = zip(*column_values)

        #columns = list(columns)

        # Get real values
        #values = [ "\"%s\"" %str(x).encode('ascii','replace').replace('\\','\\\\').replace('"','\\"') for x in values ]

        if len(special_column_values)>0:
            ncolumns, nvalues = zip(*special_column_values)
            columns.extend(ncolumns)
            values.extend(nvalues)
            del ncolumns
            del nvalues

        if on_duplicate_key is None or on_duplicate_key == "IGNORE":
            dupl_key = ""
        else:
            dupl_key = "ON DUPLICATE KEY UPDATE %s" %(on_duplicate_key)

        if on_duplicate_key == "IGNORE":
            ignore = "IGNORE"
        else:
            ignore = ""

        if self.insert_buffer is None or use_buffer==False:
            return """INSERT %s INTO %s (%s) VALUES (%s) %s""" %(ignore,
                                                                 table,
                                                                 ",".join(columns),
                                                                 ",".join(values),
                                                                 dupl_key)
        else:
            return self.insert_buffer.insert2buffer( key = self._get_buffer_key(table = table,
                                                                                columns = columns),
                                                     table = table,
                                                     columns = columns,
                                                     values = tuple(values),
                                                     max_elements_in_buffer = max_elements_in_buffer )


    def _get_buffer_multiple_queries(self, key_buffer=None):
        """
        Returns all queries of insert buffer and empties it
        
        If key_buffer is None, it generates all the queries to empty the buffer. Otherwise, it will insert only the key specified
        """
        
        return_queries = []

        if self.insert_buffer is None:
            raise ValueError("To execute insert buffer into database buffer object cannot be \"None\". First it is ncessary to create a buffer Object")
        else:
            if key_buffer is None:
                # Empties all buffer
                list_of_buffer_keys = self.insert_buffer.get_buffer().keys()
            else:
                # Empties only the key_buffer requested
                list_of_buffer_keys = [key_buffer]

            for actual_key in list_of_buffer_keys:

                bufferElement = self.insert_buffer.get_buffer()[actual_key]

                if bufferElement.num_elements>0:

                    if actual_key[0]=='I' and actual_key[1]=='U' and actual_key[2]=='_':
                        return_queries.append( self._get_multiple_insert_query( table = bufferElement.getTable(),
                                                                                columns=bufferElement.getColumns(),
                                                                                values=bufferElement.getValues(),
                                                                                externalDatabaseID = bufferElement.getSourceDBID() ) )

                    else:
                        return_queries.append( self._get_multiple_insert_query( table=bufferElement.getTable(),
                                                                                columns=bufferElement.getColumns(),
                                                                                values=bufferElement.getValues()) )

                    # It not should be here...
                    bufferElement.restart_bufferElement()


        return return_queries


    def _get_multiple_insert_query(self, table, columns, values):
        """
        "columns" must be a tuple with the name of the columns
        
        "values" must be a Set of tuples of values to insert
        """

        column_separator=", "
        values_separator=", "
        
        values_string = [ "(%s)" %(values_separator.join(x)) for x in values ]

        # tranform table objects to table name
        #table = "%s" %(table)
            
        sqlquery = "INSERT IGNORE INTO %s (%s) VALUES %s" %(str(table),
                                                            column_separator.join(columns),
                                                            column_separator.join(values_string))

        return sqlquery
    

    def _uses_buffer(self):
        """
        Returns False if it is not using buffer

        Otherwise return true
        """

        if self.insert_buffer is None:
            return False
        else:
            return True

    
    def _get_buffer_key(self, table, columns):
        """
        columns must be a list of the columns
        """
        
        return "%s%s" %(table,str(columns))


    def _lock_all_tables(self):
        
        self._unlock_tables()

        all_tables = self._get_table_names()

        self._lock_tables( table_list = all_tables )
        
        [ self.locked_tables.add(x) for x in all_tables ]
        [ self.locked_tables_win.add(x.lower()) for x in all_tables ]

        
    def _lock_tables(self, table_list=None):
        """
        Method used to Lock tables. Insertions are faster if tables are previously locked.

        table_list is a list of table names
        """

        if self.is_locked:
            raise ValueError("Trying to lock tables when there are locked tables yet")

        self.is_locked = True
        
        if table_list is None:
        	table_list = list(self.locked_tables)
        else:
        	[ self.locked_tables.add(x) for x in table_list ]
                [ self.locked_tables_win.add(x.lower()) for x in table_list ]

        if( len(table_list)>0 ):
        	self.insert_db_content( sql_query = """ LOCK TABLES %s """ %( " WRITE, ".join(table_list)+" WRITE ") )

        # COMMENTED BECAUSE VERY SLOW WHEN ADDING A NEW DATABASE WITH INSERTIONS IN LOTS OF TABLES...
        # IT SHOULD BE SUBSTITUTED BY:
                # Locking all tables
                # Put a variable telling than database is being used for insertions and then should not be used again for insertions
        #for current_autoincrement_value in self.autoincrement_values:
        #    if current_autoincrement_value[0] in self.locked_tables:
        #        self.autoincrement_values[current_autoincrement_value] = self._get_max_value(current_autoincrement_value[0],current_autoincrement_value[1])

        
    def _unlock_tables(self):
        """
        Generates the SQL statement that unlocks tables previously locked

        Method used to return the SQL query that locks access to mysql tables. Insertions are faster if tables are previously locked.

        table_list is a list of table names
        """
        
        if self._uses_buffer():
            self._empty_buffer()

        self.insert_db_content( sql_query = """UNLOCK TABLES """ )

        self.is_locked = False

    def _get_table_names(self):

        if self.table_names is None:
            self.table_names = self.select_db_content( self.GetTableNames(), answer_mode="list" )

        return self.table_names


    def _disable_indices(self, table_list=None):

        if table_list is None:
            table_list = self._get_table_names()
        else:
            if len(table_list)==0:
                table_list = self._get_table_names()

        if len(table_list)==0:
            raise ValueError("Error: there are no tables in the database.")

	[ self._check_locked_table(str(actual_table)) for actual_table in table_list ] 

        self.insert_db_content( [ "ALTER TABLE %s DISABLE KEYS" %x for x in table_list ] )
        
        return

    def _enable_indices(self, table_list=None):

        if table_list is None:
            table_list = self._get_table_names()
        else:
            if len(table_list)==0:
                table_list = self._get_table_names()

	[ self._check_locked_table(str(actual_table)) for actual_table in table_list ] 

        if len(table_list)==0:
            raise ValueError("Error: there are no tables in the database.")

        self.insert_db_content( [ "ALTER TABLE %s ENABLE KEYS" %x for x in table_list ] )
        
        return
        
    def _get_update_sql_query(self, table, update_column_values, fixed_conditions=None):
    	"""
    	"""
    	
    	if fixed_conditions is None:
    		fixed_conditions_sql = ""
    	else:
    		fixed_conditions_sql = " AND ".join( [self._get_cond_str(x) for x in fixed_conditions] )
    		
    	if fixed_conditions:
    		where_sql = " WHERE "
    	else:
    		where_sql = ""
    		
    	update_values = ",".join( ["%s=\"%s\"" %(x[0],x[1]) for x in update_column_values ] )
    		
        #print """UPDATE %s SET %s %s %s""" %(table, update_values,where_sql,fixed_conditions_sql)
    	
    	return """UPDATE %s SET %s %s %s""" %(table, update_values,where_sql,fixed_conditions_sql)
    
    	
    ####################################################################################
    #  GENERAL METHODS OF DATABASE CONTROL (LOCKS, CLOSE,...)                          #
    ####################################################################################

    def set_lock_frequency(self,frequency_value):
        """
        Method to change the lock/unlock frequency (used only in parsers, to speed up insertions and deletions
        """
        self.lock_frequency = frequency_value

        
    def _check_lock_frequency(self):

        if self.current_lock_num == self.lock_frequency:
            # Unlocking tables
            print "Unlocking and locking tables"
            self._unlock_tables()
            self.current_lock_num = 0
            self._lock_tables(self.locked_tables)
        else:
            self.current_lock_num += 1


    def _get_drop_sql_query(self, table_list):
        
        for x in table_list:
            self._check_locked_table(str(x))
        
        # Remove dropped tables from locked tables
        for x in table_list:
            self.locked_tables.discard(x)
            self.locked_tables_win.discard(x.lower())

        return ["DROP TABLE %s" %x for x in table_list ]


    def _empty_buffer(self):
        """
        It only can be used if insert buffer is being used
        """

        # First, obtain the list of all queries to insert
        queries = self._get_buffer_multiple_queries()

        for actual_query in queries:
            self.insert_db_content( actual_query, answer_mode=None )


## BUFFER RELATED METHODS ##





    
###############################################
##         BUFFER RELATED CLASSES            ##
###############################################


class Buffer(object):
    """
    Class used as an buffer for piana Inserts
    """
    def __init__(self, max_size, parent):

        self.buffer = {}

        # Maxim size for each buffer... I take a margin of 50k aprox (it works)
        #self.maxsize = int(max_size) - 50000
        self.maxsize = int(max_size) - int(int(max_size)*0.05)

        #print "Maximum buffer size: ",self.maxsize

        # Links it to the db_insert object which is using it
        self.db_insert = parent
        
        # Store the number of bytes (caracters) saved into the buffer
        #self.bytes = 0


    def get_buffer(self):
        """
        Returns the buffer dictionary
        """
        return self.buffer

    def get_buffer_element(self,key):

        if self.buffer.has_key(key):
            return self.buffer[key]
        else:
            raise ValueError("Buffer has not Buffer Element with key %s\n" %(key) )

    

    def insert2buffer( self, key, table, columns, values, max_elements_in_buffer=None ):
        """
        Inserts a query into the insert buffer

        Alert! Only used for insert buffers!!!

        Returns "None" if the query can be added to the buffer or the multiple query associated if it cannot be inserted (because buffer is full)
        """

        # Check if exists a buffer element with key "key"
        if self.buffer.has_key(key):
            # Insert if the values can be inserted in Buffer Element or restart it
            if (self.buffer[key]).insert_values(values) is None:
                multiple_query = self.db_insert._get_buffer_multiple_queries(key_buffer=key)
                self.buffer[key].restart_bufferElement(values)
                return multiple_query
        else:
            # Create a new buffer element
            self.buffer[key] = BufferElement(self.maxsize,table,columns,values,max_elements_in_buffer)

        # Control to execute one by one
        if DEBUG_BUFFER_INSERT_SINGLE:
            self.buffer[key].restart_bufferElement(values)
            return self.db_insert._get_buffer_multiple_queries(key_buffer=key)

        return None

            
class BufferElement(object):
    """
    Class used in Buffer object to store a buffer element
    """

    def __init__(self,max_size,table,columns,values=None,max_elements_in_buffer=None):
        """
        Initializes a new buffer element

        "max_size" is the maximum size of the buffer element

        "table" is the table name in the database

        "columns" is a tuple with the names of the columns

        "values" is a tuple with the values to insert

        """
        if max_size > 2000000:
            max_size = 2000000   
        self.max_size = max_size
        self.table = table
        self.columns = columns
        self.values = set()
        self.initial_size = len(table)+4

        self.size = self.initial_size

        self.max_elements_in_buffer = max_elements_in_buffer

        self.num_elements =  1  # it stores the number of elements that the buffer contains (because sometimes it is necessary to put a limit)
        

        # Sum the length of the column plus 1 (the comma)
        try:
            #print self.columns
            for i in self.columns:
                self.size += len(i)+1
        except:
            raise
            #print self.columns

        # Insert the values to the buffer
        if values:
            self.insert_values(values=values)

        #Control:
        #self.max_elements_in_buffer=1
        

    def getTable(self):

        return self.table

    def getColumns(self):

        return self.columns

    def getValues(self):

        return self.values

    def getSize(self):

        return self.size

    def insert_values(self,values):
        """
        Insert into buffer element a new tuple of values

        "values" must be a tuple of values to insert

        If the values cannot be inserted because the size has been exceeded, return None
        
        """

        # First, calculate the bytes size of this values
        # Sum the length of the value plus the comma plus the 2 parenthesis plus the 2"
        # plus the comma between two inserts
        actual_size=0

        for i in values:
            actual_size += len(str(i))+5+2

        if self.size+actual_size >= self.max_size:
            # values cannot be inserted: buffer size is exceeded
            return None
        elif self.max_elements_in_buffer is not None and self.num_elements > self.max_elements_in_buffer:
            return None
        else:
            #insert the new tuple of values and increment the size
            #self.values[values]=None
            self.values.add(values)
            self.size += actual_size

            self.num_elements += 1

            return 1


    def restart_bufferElement(self,values=None):
        """
        Restarts a Buffer Element to its initial values and insert the new tuple of values
        """
        
        # Empties the dictionary of tuples of values to insert and reset the size
                
        self.size = self.initial_size

        self.values.clear()

        self.num_elements = 0

        if values is not None:
            self.values.add(values)

            self.num_elements = 1

            # Calculate the bytes size of this values
            # Sum the length of the value plus the comma plus the 2 parenthesis plus the 2"
            # plus the comma between two inserts
            for i in values:
                self.size += len(str(i))+5+2





if __name__ == "__main__":

    print "Testing ConnectorDB"

    
    connector = DB( dbname="berit_v4", dbhost="localhost", dbuser="root", dbsocket="/home/jgarcia/local/mysql/var/mysql.sock" )
