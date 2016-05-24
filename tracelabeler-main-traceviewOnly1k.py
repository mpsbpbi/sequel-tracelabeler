import tracelabeler
import gaussmix
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
from optparse import OptionParser

parser = OptionParser("tracelabeler-windowTraceview")
parser.add_option("--trace", type="string", dest="trace")
parser.add_option("--zmw")
parser.add_option("--outputprefix")
parser.add_option("--targetFrame")
(options, args) = parser.parse_args()
if not options.trace:
    parser.print_help()
    sys.exit(1)

outprefix = options.outputprefix

tl = tracelabeler.tracelabeler( trace=options.trace )

# set the zmw and targetbase to get data
tl.setzmw(int(options.zmw))

if options.targetFrame is None:
    startframe = 0
else:
    startframe = int(options.targetFrame)

endframe = startframe+1024

#### Now traceview it to make sure everything lines up.

# plot the trace data
plt.figure(figsize=(16,2))

# green and red trace data
td = tl.traceraw(startframe,endframe)
plt.plot(range(startframe,endframe), td[0,], 'g-', linewidth=0.2)
plt.plot(range(startframe,endframe), td[1,], 'r-', linewidth=0.2)

plt.savefig('%s-traceviewOnly1k' % outprefix,dpi=500)
plt.close()

