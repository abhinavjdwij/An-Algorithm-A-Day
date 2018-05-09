'''input
'''

#~ @author = dwij28 (Abhinav Jha) ~#

class TreeNode:
	def __init__(self, x):
		self.val = x
		self.left = None
		self.right = None

class Solution:
	def inorderTraversal(self, x):
		if x == None: return []
		ans = []
		if x.left: ans += self.inorderTraversal(x.left)
		ans += [x.val]
		if x.right: ans += self.inorderTraversal(x.right)
		return ans