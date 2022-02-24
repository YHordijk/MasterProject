sigma = lambda k, n: k*n
n = 10

r = sum( sum(sigma(5, i)*sigma(3,m)*sigma(3, n-i-m) for m in range(1,n-1)) for i in range(1, n))
print(r)