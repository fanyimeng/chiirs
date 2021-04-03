def cheers(n):
    if n < 2:
        return 0
    else:
        return cheers(n - 1) + n - 1


print(cheers(1e2))

print(1e2*(1e2-1)/2)