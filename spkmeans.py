import pandas as pd
import numpy as np
import sys
import math
import warnings
import spk
warnings.simplefilter(action='ignore', category=FutureWarning)


def algorithm():
    inputs = sys.argv #inputs from cmd
    if (len(inputs) != 4):
        print("Invalid Input!")
        return # terminate
    k = int(inputs[1])
    goal = inputs[2]
    filepath = inputs[3]
    if goal != "spk": 
        info = [k, goal, filepath]
        # Use regular C method and print the output
        # saved to check if there is an error
        error = spk.standard_interface(info)
        if error == None:
            return
    else: # goal is spk
        # spk is performed using the C module
        # Get T from C:
        T = spk.spk_a([k, goal, filepath])
        if T == None: # an error
            return
        df = np.array(T)
        k = len(df[0])
        d = k
        keys = create_cents(df, k, d)  # perform kmeans++
        info = [keys, T]
        # print the selected keys:
        for value in keys[:-1]:
            print(value, end=",")
        print(keys[-1])
        # Run second part of spk
        spk.spk_b(info)


def create_cents(df, k, d):
    N = len(df)
    np.random.seed(0)
    indexes = [i for i in range(N)]
    # get the first value rendomaly
    index1 = np.random.choice(indexes)
    # list of keys
    keys = [index1]
    # first centroid
    mu1 = df[index1]
    cents = [mu1]
    # kmeans++ algo:
    Z = 1
    while Z < k:
        d_s = []
        for i in range(N):
            d_i = math.inf
            x_i = df[i]
            for j in range(Z):
                arg_min = 0
                for p in range(d):
                    arg_min = arg_min+(float(x_i[p])-float(cents[j][p]))**2
                if arg_min < d_i:
                    # replace the minimum
                    d_i = arg_min
            # add d_i
            d_s.append(d_i)
        # calculate probs
        probs = []
        for i in range(N):
            probs.append(d_s[i]/sum(d_s))
        index_z = np.random.choice(indexes, p=probs)
        keys.append(index_z)
        mu_z = df[index_z]
        cents.append(mu_z)
        Z += 1

    results = []
    for cent in cents:
        results.append(cent)
    return keys


if __name__ == "__main__":
    algorithm()
