## test SCRIPT
import datetime

import sys
sys.path.insert(1, './py_functions/')
from time_functions import datenum2date, date2datenum


datenum=719529.5
d=datetime.datetime(1970,1,1,12)

print(datenum2date(datenum))
print(date2datenum(d))
