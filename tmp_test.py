## test SCRIPT
import datetime
import sys
sys.path.insert(1, './py_functions/')
from time_functions import datenum2date, date2datenum


datenum=719529.5374
d=datetime.datetime(1970,1,1,12,3,4,3)

dt2d = datenum2date(datenum)
print(dt2d)
d2dt = date2datenum(d)
print(d2dt)
