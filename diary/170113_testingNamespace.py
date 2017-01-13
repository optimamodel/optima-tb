
try:
    import optima
    print("Could import optima")
except:
    print("couldn't import optima")
    
try:
    import optima.tb
    print("Could import optima.tb")
except:
    print("couldn't import optima.tb")
    
   
try:
    import optima.tb.utils
    print("Could import optima.tb.utils")
except:
    print("couldn't import optima.tb.utils")

try:
    from optima.tb import utils
    print("Could 'from optima.tb import utils'")
except:
    print("couldn't 'from optima.tb import utils'")
        
   
try:
    import optima.tb.utils as utils
    print("Could import optima.tb.utils as utils")
except:
    print("couldn't import optima.tb.utils as utils")
    
    
from optima.tb import settings    
import optima.tb.settings
from optima.tb.project import Project
import optima.tb.settings as settings
