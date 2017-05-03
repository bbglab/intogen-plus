from bgparsers import readers, selector
from pprint import pprint

pancanatlas = '/home/jordeu/workspace/intogen/intogen-plus/pancanatlas'

for mut in readers.variants(selector.find(pancanatlas)):
    pprint(mut)
    break
print("\n")