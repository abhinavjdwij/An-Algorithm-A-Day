'''input
2 6
'''

#~ @author = dwij28 (Abhinav Jha) ~#

def expo(a, b):
    if a == 0: return 0
    elif b == 0: return 1
    elif b % 2 == 0: return expo(a, b/2) * expo(a, b/2)
    return a * expo(a, b/2) * expo(a, b/2)

number, power = map(int, raw_input().rstrip().split())
print expo(number, power)