#!/usr/bin/env python3

import sys

input = open(sys.argv[1], 'r')
start_city = input.readline().rstrip()
cities = set()
for i, l in enumerate(input):
    [f, t, d, p] = l.rstrip().split(' ')
    cities.add(f)
input.close()

print("{} {}".format(len(cities), start_city))
for i, city in enumerate(cities):
    print(i)
    print(city)
input = open(sys.argv[1], 'r')
input.readline()
for l in input:
    [f, t, d, p] = l.rstrip().split(' ')
    print("{} {} {} {}".format(f, t, int(d)+1, p))
input.close()
