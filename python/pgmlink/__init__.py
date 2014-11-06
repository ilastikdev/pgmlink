import pgmlink as _pgmlink
from TransitionClassifier import *

#Load all variables of pgmlink to the top level
#so as to allow accessing them by 'pgmlink.Traxel'
for i in dir(_pgmlink):
    if i.startswith("__"):
        continue
    exec(i+" =  _pgmlink."+i)