'''input
5
1 2 3 3 5
'''

#~ @author = dwij28 (Abhinav Jha) ~#

def solve(arr):
    s, ssq, n = 0, 0, len(arr)
    for i in arr:
        s, ssq = s+i, ssq+(i*i)
    a_minus_b = s - (n*(n+1)/2)
    asq_minus_bsq = ssq - (n*(n+1)*(n+n+1)/6)
    a_plus_b = asq_minus_bsq / a_minus_b
    a = (a_plus_b + a_minus_b) / 2
    b = (a_plus_b - a_minus_b) / 2
    return [a, b]

n, data = input(), map(int, raw_input().strip().split())
repeat, missing = solve(data)
print repeat, 'appears twice and', missing, 'does not appear.'