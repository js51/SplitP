# Norm manually
start = time.time()
rows, cols = flattening.nonzero()
f_norm_1 = sum(flattening[i,j]**2 for i, j in zip(rows, cols))**(1/2)
end = time.time()
print(end-start)
# Builtin
start = time.time()
f_norm_2 = norm(flattening, 'fro')
end = time.time()
print(end-start)
# Compare
print(f_norm_1 == f_norm_2) # True
