from SphinxReport.Tracker import *

import CGAT.Pipeline as P

#############################################################################
# Get parameterization

P.getParameters( 
    ["%s/pipeline.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )

PARAMS = P.PARAMS
class ProjectTracker(TrackerSQL):


    PARAMS = P.PARAMS

    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = "sqlite:///./"+self.PARAMS['database'] , **kwargs )
        
        # issuing the ATTACH DATABASE into the sqlalchemy ORM (self.db.execute( ... ))
        # does not work. The database is attached, but tables are not accessible in later
        # SELECT statements.
        if not self.db:
            def _create():
                import sqlite3
                conn = sqlite3.connect(re.sub( "sqlite:///", "", self.PARAMS['database']) )
                statement = '''ATTACH DATABASE '%s' as annotations;
                   ATTACH DATABASE 'mapping.dir/csvdb' as mapping; ''' % self.PARAMS["annotations_database"]
                                                             
         
     
                conn.executescript(statement)
                print statement 
                return conn
              
            self.connect( creator = _create )
