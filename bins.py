import os, sys

bin_count = 100
if len(sys.argv) > 1:
  bin_count = int(sys.argv[1])

bins = [[] for _ in xrange(bin_count)]

with open('distances2.txt') as f:
  for line in f:
    num = float(line)
    idx = int(num * bin_count / .1)
    if idx >= bin_count: continue
    bins[idx].append(num)

for idx, b in enumerate(bins):
  print "%.9f --> %d %s" % (idx * .1 / bin_count, len(b), b)
