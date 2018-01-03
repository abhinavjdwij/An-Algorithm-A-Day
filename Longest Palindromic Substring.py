'''input
abbabcdabbabbbadb
'''

#~ @author = dwij28 (Abhinav Jha) ~#

'''

************************ Resources ************************

1: Brute Force - O(n^3) and Dynamic Programming - O(n^2) Time + O(n^2) Space Approach: https://www.geeksforgeeks.org/?p=25463
2. O(n^2) Time and O(1) Space Approach - O(n^2): https://www.geeksforgeeks.org/longest-palindromic-substring-set-2/
3. Manacher's Algorithm - O(n): https://www.youtube.com/watch?v=V-sEwsca1ak

'''

def lps(s): # Length of Longest Palindromic Substring
	ans, n, l, r = 1, len(s), 0, 0
	for i in xrange(1, n):
		l, r = i-1, i
		while l >= 0 and r < n and s[l] == s[r]:
			ans = max(ans, r-l+1)
			l, r = l-1, r+1
		l, r = i-1, i+1
		while l >= 0 and r < n and s[l] == s[r]:
			ans = max(ans, r-l+1)
			l, r = l-1, r+1
	return ans

	# O(n^2)

s = raw_input().strip()
print lps(s)