
OPTIMIZE_SPACE = None


# FILE TO IMPLEMENT THE CLASSES RELATED WITH A DATABASE
class SQLSelectStatement(object):

    def __init__(self):

        self.tables = set()
        self.columns = []
        self.fixed_conditions = set()
        self.join_conditions = set()
        self.group_conditions = set()



    def add_element(self, columns=[], tables=[], fixed_conditions = [], join_conditions = [], group_conditions = [] ):

        [ self.tables.add(x) for x in tables ]
        [ self.columns.append(x) for x in columns ]
        [ self.fixed_conditions.add(x) for x in fixed_conditions ]
        [ self.join_conditions.add(x) for x in join_conditions ]
        [ self.group_conditions.add(x) for x in group_conditions ]

    def reset_columns(self, columns = []):

        self.columns = columns

    def merge(self,sqlSelectStat):

        self.tables.update(sqlSelectStat.tables)
        self.columns.extend(sqlSelectStat.columns)
        self.fixed_conditions.update(sqlSelectStat.fixed_conditions)
        self.join_conditions.update(sqlSelectStat.join_conditions)
        self.group_conditions.update(sqlSelectStat.group_conditions)

    def __repr__(self):
	return "Tables: %s\nColumns: %s\nFixed Conditions: %s\nJoin Conditions: %s\nGroup Conditions: %s\n" % (map(str, self.tables), self.columns, self.fixed_conditions, self.join_conditions, self.group_conditions)


class Database(object):
    """
    Class to represent a relational database
    """

    def __init__(self, tables=None):
        """
        "tables": a list of table objects that are contained in the database
        """
        
        if tables is None:
            tables = []

        self.tables = tables
        self.created_unique_tables = {}  # stores the unique tables created ( the unique_table_name as a key and the optimization number as value )
        self.added_unique_table_names = set()  # Stores the unique tables that have been previously inserted in the tables list
        self.all_tables = set()          # stores the name of all tables, unique and temporal tables included

    def remove_tables(self):
        self.tables = []
        self.created_unique_tables = {}  # stores the unique tables created ( the unique_table_name as a key and the optimization number as value )
        self.added_unique_table_names = set()  # Stores the unique tables that have been previously inserted in the tables list
        self.all_tables = set()

    def get_sql_query( self, ignore_primary_keys = False ):
        """
        """
        sql_query_string = "\n".join( [x.create_mysql_query(self.created_unique_tables, ignore_primary_key=ignore_primary_keys ) for x in self.tables] )
        self.created_unique_tables = {}
        return sql_query_string


    def add_table(self, table):

        if table.get_table_name() in self.all_tables:
            raise ValueError("Trying to insert the same table twice")

        self.tables.append(table)
        self.all_tables.add(table.get_table_name())

        for actual_field in table.get_fields():
            if actual_field.get_optimized_space():
                if actual_field.get_optimized_table().get_table_name() not in self.added_unique_table_names:
                    self.all_tables.add(actual_field.get_optimized_table())
                    self.added_unique_table_names.add(actual_field.get_optimized_table().get_table_name())
                self.all_tables.add(table.get_temp_table())


    def get_tables(self):
        """
        Get basic tables
        """
        #return list(self.all_tables)
        return self.tables

    def get_all_tables(self):

        return list(self.all_tables)

    def get_drop_sql_query(self):

        return "\n".join([x.get_drop_query() for x in self.tables])


    def optimize_database(self, optimize=False, analyze=False):
        """
        Returns the queries to optimize database
        """

        queries = set()

        if analyze:
            for actual_table in self.tables:
                queries.update(actual_table.get_analyze_query())


        if optimize:
            for actual_table in self.tables:
                queries.update(actual_table.get_optimize_queries())

        return list(queries)


class FieldDB(object):
    """
    Class to represent a field in a database
    """

    def __init__(self, field_name, data_type, default_value=None, null=True, foreign_key=None, user_friendly_name=None, optimize_space = None, compress = None ):
        """
        "field_name": field name in the database
        "data_type": data type of this field
        "default_value"
        "null": By default, it is accecpted a null value. If not, it is necessary to specify it
        "foreign_key": it is a tuples with the table and the Field of its correspondent key
                        (TableObject, FieldObject)
                        They should be of the same data_type!!!
        "user_friendly_name": This is used to refer to this field in a user friendly manner (for different fields in different tables that have the same meaning...)
        "optimize_space": This is used to apply the space optimization for this field. It consists on saving the unique values of this field in a separated table, and an integer cross-reference is saved in the original table. It can be None, or the number of bytes reserved to index it
        """

        self.field_name = field_name
        self.data_type = data_type
        self.default_value = default_value
        self.null = null
        self.foreign_key = foreign_key

        if OPTIMIZE_SPACE:
            self.optimize_space = optimize_space
        else:
            self.optimize_space = None

        self.optimized_table = None

        self.compress = compress

        if self.optimize_space is not None:
            self.optimized_table = self.get_optimized_table()

        if( user_friendly_name is None ):
            self.user_friendly_name = None
        else:
            self.user_friendly_name = user_friendly_name.lower()

        self.sqlSelect = None

    def is_compressed(self):
        return self.compress

    def get_optimized_field_name(self):
        return "original_id_%s" %(self.field_name)
        
    def get_optimized_space(self):
        """
        Returns the number of bytes for the optimized form of the field, if it has. Otherwise, return None
        """
        return self.optimize_space
        
    def get_user_friendly_name(self):

        if self.user_friendly_name is None:
            return self.field_name

        return self.user_friendly_name

    def get_field_name(self):
        """

        """
        return self.field_name

    def get_data_type(self):
        """

        """
        return self.data_type

    def get_optimized_data_type(self):
        return "integer(%s) unsigned" %(self.optimize_space)

    def get_mysql_query(self):

        if self.optimize_space is None:
            data_type = self.get_data_type()
        else:
            data_type = self.get_optimized_data_type()
        
        query = "%s %s " %(self.get_field_name(),
                           data_type)

        if self.default_value is not None:
            query += " DEFAULT \"%s\"" %(self.default_value)
            
        if self.null is False:
            query += " NOT NULL"

        return query

    # It can be done more elegant...
    def get_not_optimized_query(self):

        if self.optimize_space is None:
            return self.get_mysql_query()

        query = "%s %s " %(self.get_optimized_field_name(),
                           self.get_data_type())

        if self.default_value is not None:
            query += " DEFAULT \"%s\"" %(self.default_value)
            
        if self.null is False:
            query += " NOT NULL"

        return query


    def set_foreign_key(self, foreign_key):
	"""
	"""

        self.foreign_key = foreign_key


    def get_optimized_table(self):
        """
        Get the optimized table
        """

        if self.optimize_space is None:
            raise ValueError("Trying to get an optimized table in a not allowed field")

        if self.optimized_table is None:

            self.optimized_table = TableDB( table_name = "unique_%s" %(self.field_name),
                                            table_fields = [ FieldDB( field_name = self.get_field_name(),
                                                                      data_type = "%s auto_increment" %(self.get_optimized_data_type()),
                                                                      null = True),
                                                             FieldDB( field_name = self.get_optimized_field_name(),
                                                                      data_type = self.get_data_type(),
                                                                      null = True) ],
                                            primary_key = (self.get_optimized_field_name()),
                                            indices = [(self.get_field_name())]  )

        return self.optimized_table

##     def get_optimized_temp_table(self):
##         """
##         Get the temp table to be able to generate quickly the optimized data
##         """

##         if self.optimize_space is None:
##             raise ValueError("Trying to get an optimized temp table in a not allowed field")

##         if self.optimized_temp_table is not None:
##             return self.optimized_temp_table
##         else:

##             new_table = TableDB( table_name = "temp_%s" %(self.field_name),
##                                  table_fields = [ FieldDB( field_name = self.get_field_name(),
##                                                            data_type = "%s auto_increment" %(self.get_optimized_data_type()),
##                                                            null = 1),
##                                                   FieldDB( field_name = self.get_optimized_field_name(),
##                                                            data_type = self.get_data_type(),
##                                                            null = 1) ],
##                                  primary_key = (self.get_optimized_field_name()),
##                                  indices = [(self.get_field_name())]  )

##             return new_table
        


    def get_foreign_key_mysql_query(self):

        if self.foreign_key is not None:
            return "FOREIGN KEY (%s) REFERENCES %s(%s)" %(self.field_name,
                                                          self.foreign_key[0].get_table_name(),
                                                          self.foreign_key[1].get_field_name())
                          


class TableDB(object):
    """
    This object will implement a database table
    """

    def __init__(self, table_name, table_fields, primary_key=(), indices=[], fulltext_indices = []):
        """
        Intializes the object

        "table_name" will be the table name in the database
        "table_fields" will be a list of fieldDB Objects
        "primary_key" A tuple or list of the name of the fields that form the primary_key
        "indices": a list of tuples with the name of the fields that form the indices
        """

        self.table_name = table_name
        self.table_fields = table_fields
        self.primary_key = primary_key
        self.indices = indices
        self.fulltext_indices = fulltext_indices

        self.temp_table_name = None
        self.has_optimized_field = self._has_optimized_fields()
        self.temp_table = None


        # Stores the information of fields by field_name
        self.fields_by_name = {}
        for actual_field in self.table_fields:
            self.fields_by_name[actual_field.get_field_name()] = actual_field

        self.sqlSelectStats = {}

    def __str__(self):
        return self.table_name

    def set_primary_key(self, primary_key):
        self.primary_key = primary_key

    def has_optimized_fields(self):
        return self.has_optimized_field

    def _has_optimized_fields(self):

        for x in self.table_fields:
            if x.get_optimized_space() is not None:
                self.temp_table_name = "temp_"+self.get_table_name()
                return 1

        return None
            
    def get_temp_table_name(self):
        return self.temp_table_name

    def add_field(self, new_field):
        self.table_fields.append(new_field)

    def has_field(self, field_name):

        if self.fields_by_name.has_key(field_name):
            return 1
        else:
            return None


    def get_fields(self):

        return self.table_fields
        

    def get_table_name(self):

        return self.table_name

    def set_table_name(self, new_name):

        self.table_name = new_name


    def get_field(self, field_name):

        if self.fields_by_name.has_key(field_name):
            return self.fields_by_name[field_name]
        else:
            return None

    def add_indice(self, indice):

        self.indices.append(indice)


    def has_indice(self, indice):
        """
        Checks if the table contains an indice with those fields
        """
        if( self.indices.contains(indice) ):
            return 1
        else:
            return None


    def get_optimize_queries(self):

        query_list = [ "OPTIMIZE TABLE %s" %self.table_name ]
        
        for x in self.table_fields:
            if x.get_optimized_space() is not None:
                query_list.append( "OPTIMIZE TABLE %s" %x.get_optimized_table().get_table_name() )
                
        return query_list

    def get_analyze_query(self):

        query_list = [ "ANALYZE TABLE %s" %self.table_name ]

        for x in self.table_fields:
            if x.get_optimized_space() is not None:
                query_list.append( "ANALYZE TABLE %s" %x.get_optimized_table().get_table_name() )

        return query_list

    def get_SQLSearchStatement(self, field, value_list):

        sqlStat = SQLSelectStatement()
        sqlStat.add_element(tables = [self.table_name])

        field = self.get_field(field)
        if field is None:
            raise ValueError("Trying to get a field (%s) that is not found in thid table" %field)
        if field.get_optimized_space() is not None:
            sqlStat.add_element(tables=[field.get_optimized_table().get_table_name()])
            sqlStat.add_element(join_conditions = [("%s.%s" %(self.table_name,field.get_field_name()),
                                                    "=",
                                                    "%s.%s" %(field.get_optimized_table().get_table_name(),field.get_field_name()) ) ])
            #sqlStat.add_element(fixed_conditions = [(field.get_optimized_field_name(),"=",value)])
            sqlStat.add_element(fixed_conditions = [(field.get_optimized_field_name(),"IN","(\"%s\")" %("\",\"".join([ str(x) for x in value_list ])),None)])
        else:
            sqlStat.add_element(fixed_conditions = [(field.get_field_name(),"IN","(\"%s\")" %("\",\"".join([ str(x) for x in value_list ])),None)])
            #sqlStat.add_element(fixed_conditions = [(field.get_field_name(),"=",value)])

        return sqlStat

    def get_SQLSelectStatement(self, field_list):

        if not self.sqlSelectStats.has_key(tuple(field_list)):

            sqlStat = SQLSelectStatement()
            sqlStat.add_element(tables=[self.table_name])

            for actual_field in field_list:
                field = self.get_field(actual_field)
                if field is None:
                    raise ValueError("Trying to get a field (%s) that is not found in this table" %(actual_field))
                if field.get_optimized_space() is not None:
                    sqlStat.add_element(tables=[field.get_optimized_table().get_table_name()])
                    sqlStat.add_element(join_conditions = [("%s.%s" %(self.table_name,field.get_field_name()),
                                                            "=",
                                                            "%s.%s" %(field.get_optimized_table().get_table_name(),field.get_field_name()) ) ])
                    sqlStat.add_element(columns = ["%s.%s" %(field.get_optimized_table().get_table_name(),field.get_optimized_field_name())])
                else:
                    sqlStat.add_element(columns = ["%s.%s" %(self.get_table_name(),field.get_field_name())])

            self.sqlSelectStats[tuple(field_list)] = sqlStat

        return self.sqlSelectStats[tuple(field_list)]

    
    def create_mysql_query(self, created_unique_tables={}, ignore_primary_key=False):

        sql_fields = [ x.get_mysql_query() for x in self.table_fields ]

        
        if ignore_primary_key is False:
            # In this case, the primary key is saved as a normal index
            primary_key_text = "PRIMARY KEY"
        else:
            primary_key_text = "INDEX"


        if isinstance(self.primary_key,list) or isinstance(self.primary_key,tuple):
            if len(self.primary_key)>0:
                primary_key_sql = "%s (%s)" %(primary_key_text, ",".join(list(self.primary_key)))
            else:
                primary_key_sql = None
        else:
            primary_key_sql = "%s (%s)" %(primary_key_text, self.primary_key)

        if primary_key_sql is not None:
            sql_fields.append(primary_key_sql)

        indices_list = []
        for actual_indice in self.indices:
            if( isinstance(actual_indice,tuple) ):
                indices_list.append([ x for x in actual_indice ])
            else:
                indices_list.append([actual_indice])
               
        indices_sql = ["INDEX (%s)" %(",".join(actual_indice)) for actual_indice in indices_list]

        sql_fields.extend(indices_sql)

        fulltext_list = []
        for current_index in self.fulltext_indices:
            if( isinstance(current_index,tuple) ):
                fulltext_list.append([ x for x in current_index ])
            else:
                fulltext_list.append([current_index])

        sql_fields.extend(["FULLTEXT (%s)" %",".join(current_fulltext_index) for current_fulltext_index in fulltext_list ])


        # ADD FOREIGN KEY RESTRICTIONS
        for actual_field in self.table_fields:
            foreign_key_query = actual_field.get_foreign_key_mysql_query()
            if foreign_key_query is not None:
                sql_fields.append(actual_field.get_foreign_key_mysql_query())

        engine = "MyISAM"

        query = ["CREATE TABLE IF NOT EXISTS %s ( \n\t%s\n ) ENGINE %s;" %(self.table_name,
                                                             ",\n\t".join(sql_fields),
                                                             engine)]

        # CREATE  AUXILIAR TABLES

        if self.has_optimized_field is not None:
            #indices = []
            for actual_field in self.table_fields:
                if actual_field.get_optimized_space() is not None:
                    opt_table = actual_field.get_optimized_table()
                    if not created_unique_tables.has_key(opt_table.get_table_name()):
                        query.append(actual_field.get_optimized_table().create_mysql_query(created_unique_tables))
                        created_unique_tables[opt_table.get_table_name()] = 1
                    #indices.append(actual_field.get_optimized_field_name())
        
            sql_fields = [ x.get_not_optimized_query() for x in self.table_fields ]

            #indices_sql = ["INDEX (%s)" %actual_indice for actual_indice in indices]
            #sql_fields.extend(indices_sql)
                                      
            
            query.append("CREATE TABLE %s ( \n\t%s\n ) ENGINE MyISAM;" %(self.get_temp_table_name(),
                                                                         ",\n\t".join(sql_fields)) )


        return "\n".join(query)


    def get_temp_table(self):
        """
        """
        
        if self.has_optimized_field is None:
            return None

        if self.temp_table is None:
            self.temp_table = TableDB( table_name = self.get_temp_table_name(),
                                       table_fields = self.get_fields() )
        
        return self.temp_table


    def get_drop_query(self):

        return "DROP TABLE %s;" %self.table_name


    def empty_temporal_table_sql_query(self, _skip_first=None):

        if self.has_optimized_field is None:
            return None


        # Steps:

        queries = []

        # 1. Put all distinct values into its unique tables
        if _skip_first is None:
            for actual_field in self.table_fields:
                if actual_field.get_optimized_space() is not None:
                    queries.append("ALTER TABLE %s DISABLE KEYS" %(actual_field.get_optimized_table()) )
                    queries.append("INSERT IGNORE INTO %s (%s) SELECT DISTINCT %s FROM %s" %(actual_field.get_optimized_table(),
                                                                                    actual_field.get_optimized_field_name(),
                                                                                    actual_field.get_optimized_field_name(),
                                                                                    self.get_temp_table_name() ) )
                    queries.append("ALTER TABLE %s ENABLE KEYS" %(actual_field.get_optimized_table()) )


        # 2. Transfer all the information from the temp table to the final table
        unique_tables = []
        values_in_temp_table = []
        where_conditions = []
        for actual_field in self.table_fields:
            if actual_field.get_optimized_space() is not None:
                unique_tables.append( actual_field.get_optimized_table().get_table_name() )
                values_in_temp_table.append( "%s.%s" %(actual_field.get_optimized_table(),
                                                       actual_field.get_field_name()))
                where_conditions.append( "%s.%s = %s.%s" %(actual_field.get_optimized_table(),actual_field.get_optimized_field_name(),
                                                           self.get_temp_table_name(),actual_field.get_optimized_field_name()))
            else:
                values_in_temp_table.append( actual_field.get_field_name() )
                                                       
        queries.append("ALTER TABLE %s DISABLE KEYS" %(self.get_table_name()))
        queries.append( "INSERT IGNORE INTO %s (%s) SELECT %s FROM %s, %s WHERE %s" %(self.get_table_name(),
                                                                                      ",".join( [x.get_field_name() for x in self.table_fields] ),
                                                                                      ",".join(values_in_temp_table),
                                                                                      self.get_temp_table_name(),
                                                                                      ",".join(unique_tables),
                                                                                      " AND ".join(where_conditions) ) )
        queries.append("ALTER TABLE %s ENABLE KEYS" %(self.get_table_name()))
        
        # 3. Empty the temp table
        queries.append( "DELETE FROM %s" %self.get_temp_table_name() )

        #print queries
        return queries
        
    def __repr__(self):
	return "%s" % self.table_name

