'''input
daily algorithm practice
'''

#~ @author = dwij28 (Abhinav Jha) ~#

def solve(s):
	data = s.split()[::-1]
	ans = ' '.join(data)
	return ans

s = raw_input().strip()
print solve(s)