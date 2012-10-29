import sys
import dspdsmks

params = dspdsmks.Parameters(sys.argv[1])
for p in params:
    print "%s & %s & %s\\" % (p.name, p.value, p.error)
