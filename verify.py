#!/usr/bin/env python3

import numpy as np
import os
import subprocess
import sys
import time


def verify(input_filename, output_filename, time = None, silent=False):
    try:
        input = open(input_filename, 'r')
        output = open(output_filename, 'r')
    except:
        print("Could not open files")
        sys.exit(1)

    input_cities = set()
    input_flights = dict()
    input_start_city = ''

    regions = []
    input_region_count, input_start_city = input.readline().rstrip().split(" ")
    input_region_count = int(input_region_count)
    for i in range(input_region_count):
        region_name = input.readline()
        cities = input.readline().rstrip().split(" ")
        regions.append(cities)
        if input_start_city in cities:
            input_start_region = cities
        for c in cities:
            input_cities.add(sys.intern(c))

    for l in input:
        [f, t, d, p] = l.rstrip().split(' ') #From_city, To_city, Day, Price
        f = sys.intern(f)
        t = sys.intern(t)
        d = int(d)
        p = int(p)

        days = [d] if d > 0 else range(1, input_region_count+1)
        for d in days:
            if (f, t, d) in input_flights:
                if p < input_flights[(f, t, d)]:
                    input_flights[(f, t, d)] = p
            else:
                input_flights[(f, t, d)] = p
        input_cities.add(f)
        input_cities.add(t)
    input.close()

    output_flights = list()
    output_price = 0
    calculated_price = 0

    try:
        output_price = int(output.readline().rstrip())
        for l in output:
            [f, t, d, p] = l.rstrip().split(' ') #From_city, To_city, Day, Price
            f = sys.intern(f)
            t = sys.intern(t)
            d = int(d)
            p = int(p)

            if (f, t, d) not in input_flights:
                if p < 500000 or not silent:
                    print("Error: Flight on output is not in input: {} {} {} {}".format(f, t, d, p))
                return (999999, 0)
            calculated_price += input_flights[(f, t, d)]
            output_flights.append((f, t, d))
        output.close()
    except:
        print("Error: Wrong input -- error on the very first line")
        return (999999, 0)

    if calculated_price != output_price:
        print("Error: Price is not good: is {}, should be {}".format(output_price, calculated_price))

    # check if flight sequence is good
    for i in range(len(output_flights)-1):
        if output_flights[i][1] != output_flights[i+1][0]:
            print("Error: Sequence of flights is not good")

    #check if starting city is good
    if input_start_city != output_flights[0][0]:
        print("Error: Starting city is wrong")

    #check ending city
    if output_flights[-1][1] not in input_start_region:
        print("Error: Ending city is not in starting region ({})".format(", ".join(input_start_region)))

    #check that each flight departures on good day
    for i in range(len(output_flights)):
        if output_flights[i][2] != i+1:
            print("Flights departures on wrong day")

    #check that each city is visited
    if len(output_flights) != input_region_count:
        print("Not each region is visited")

    standards = {
        "TSALESMAN2-1.in":   1396,
        "TSALESMAN2-2.in":   2159,
        "TSALESMAN2-3.in":  44151,
        "TSALESMAN2-4.in": 118814,
        
        "data_5.txt.in":     1950,
        "data_10.txt.in":    5375,
        "data_15.txt.in":    4281,
        "data_20.txt.in":    6053,
        "data_30.txt.in":    7629,
        "data_40.txt.in":    7660,
        "data_50.txt.in":    7308, # (ve fóru 7235)
        "data_60.txt.in":    9099, # (ve fóru 9180)
        "data_70.txt.in":   11907, # (ve fóru 12358)
        "data_100.txt.in":  14140, # (ve fóru 14942)
        "data_200.txt.in":  28124, # (ve fóru 28338, ~27000)
        "data_300.txt.in":  41235, # (ve fóru 34517)

        "pidi.in":              4,
        "mini.in":            100,
        "micro.in":           100,
        "pico.in":            101,
        "mili.in":            111,
        "nano.in":             60,
    }
    standard_price = 0
    for k, v in standards.items():
        if input_filename.endswith(k):
            standard_price = v
            break
    weight = 1.0
    if input_filename.endswith("TSALESMAN2-3.in") or input_filename.endswith("TSALESMAN2-4.in"):
        weight = 1.5

    if not silent:
        print("OK -- price: {}, points: {:.2f}, weighted: {:.2f}{}".format(
            output_price,
            standard_price * 100.0 / output_price,
            standard_price * 100.0 / output_price * weight,
            "" if not time else ", time: {:.3f}".format(time),
        ))
    return (output_price, standard_price * 100.0 / output_price * weight)


def test(source_filename, input_filename, runs):
    print("Compiling ... ", end="")
    sys.stdout.flush()
    os.system("./compile.sh {}".format(source_filename))
    print("OK")
    filename, extension = os.path.splitext(os.path.split(input_filename)[1])
    output_filename = os.path.join("output", filename + ".out")

    times = []
    scores = []
    prices = []
    print("Runninng and verifying ... ", end="")
    for i in range(int(runs)):
        print(i+1, end=" ")
        sys.stdout.flush()
        time_start = time.time()
        os.system("./a.out < {} > {}".format(input_filename, output_filename))
        time_end = time.time()
        price, score = verify(input_filename, output_filename, silent=True)
        prices.append(price)
        scores.append(score)
        times.append(time_end - time_start)
    print("OK")
    print("Price: mean {:10.3f}, min {:7d},     max {:7d}".format(np.mean(prices), int(np.min(prices)), int(np.max(prices))))
    print("Score: mean {:10.3f}, min  {:10.3f}, max  {:10.3f}".format(np.mean(scores), np.min(scores), np.max(scores)))
    print("Time:  mean {:10.3f}, min      {:6.3f}, max      {:6.3f}".format(np.mean(times), np.min(times), np.max(times)))

if len(sys.argv) == 2:
    directory = sys.argv[1]
    score = 0
    for f in sorted(os.listdir(directory)):
        if not f.endswith(".in") and not f.endswith(".dat"):
            continue
        filename, extension = os.path.splitext(f)
        print("=== {} ===".format(filename))
        fin = os.path.join(directory, f)
        fout = os.path.join("output", filename + ".out")
        time_start = time.time()
        os.system("./a.out < {} > {}".format(fin, fout))
        time_end = time.time()
        score += verify(fin, fout, time_end - time_start)[1]
    print("\nSum: {:.2f}".format(score))
elif len(sys.argv) == 3:
    verify(sys.argv[1], sys.argv[2])
elif len(sys.argv) == 4:
    test(sys.argv[1], sys.argv[2], sys.argv[3])
else:
    print("""Expecting
either 1 argument - 1) Input directory,
or 2 arguments - 1) Input file, 2) Contestant\'s output file
or 3 arguments - 1) Source code file, 2) Input file, 3) Number of runs
""")
    sys.exit(1)
