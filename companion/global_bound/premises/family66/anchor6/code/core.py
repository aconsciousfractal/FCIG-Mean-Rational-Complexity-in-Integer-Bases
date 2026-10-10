"""Only the fresh marked-6 anchor is supported."""
from exact import need
def factors(p):
 need(p==1,"marked-6 sentinel");return (2,3)
def mark(p):
 factors(p);return 6,12
