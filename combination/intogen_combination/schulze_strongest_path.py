def strongest_path(size, pref, spath):

    for i in range(size):
        for j in range(size):
            if i != j:
                if pref[i*size + j] > pref[j*size + i]:
                    spath[i*size + j] = pref[i*size + j]

    for i in range(size):
        for j in range(size):
            if i != j:
                for k in range(size):
                    if (i != k) and (j != k):
                        spath[j*size + k] = max(spath[j*size + k], min(spath[j*size + i], spath[i*size + k]))