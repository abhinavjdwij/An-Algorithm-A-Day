'''input
A+(B*C-(D/E-F)*G)*H
'''

#~ @author = dwij28 (Abhinav Jha) ~#

op = {'+' : 1, '-' : 1, '*' : 2, '/' : 2, '^' : 3}

def prec(x): return 0 if x not in op else op[x]

def infixToPrefix(s):
	s = list(s[::-1])
	for i in xrange(len(s)):
		if s[i] == '(': s[i] = ')'
		elif s[i] == ')': s[i] = '('
	ans, stack = [], []
	for i in s:
		if i in op:
			while stack and prec(i) <= prec(stack[-1]):
				ans.append(stack[-1])
				stack.pop()
			stack.append(i)
		elif i == '(': stack.append(i)
		elif i == ')':
			while stack and stack[-1] != '(':
				ans.append(stack[-1])
				stack.pop()
			if stack and stack[-1] == '(': stack.pop()
		else: ans.append(i)
	while stack:
		if stack[-1] not in '()': ans.append(stack[-1])
		stack.pop()
	ans = ans[::-1]
	return ''.join(ans)

s = raw_input().strip()
print infixToPrefix(s)