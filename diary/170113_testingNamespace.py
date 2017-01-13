
try:
    import optima
    print("Could import optima")
except:
    print("couldn't import optima")
    
try:
    import optima_tb
    print("Could import optima_tb")
except:
    print("couldn't import optima_tb")
    
   
try:
    import optima_tb.utils
    print("Could import optima_tb.utils")
except:
    print("couldn't import optima_tb.utils")

try:
    from optima_tb import utils
    print("Could 'from optima_tb import utils'")
except:
    print("couldn't 'from optima_tb import utils'")
        
   
try:
    import optima_tb.utils as utils
    print("Could import optima_tb.utils as utils")
except:
    print("couldn't import optima_tb.utils as utils")
    
    
from optima_tb import settings    
import optima_tb.settings
from optima_tb.project import Project
import optima_tb.settings as settings
