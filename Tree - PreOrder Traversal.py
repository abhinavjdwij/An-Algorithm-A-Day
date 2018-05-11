'''input
'''

#~ @author = dwij28 (Abhinav Jha) ~#

class TreeNode:
	def __init__(self, x):
		self.val = x
		self.left = None
		self.right = None

def postorderTraversal(self, x):
    ans = []
    if x == None: return []
    if x.left: ans += self.postorderTraversal(x.left)
    if x.right: ans += self.postorderTraversal(x.right)
    ans += [x.val]
    return ans