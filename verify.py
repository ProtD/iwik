#!/usr/bin/env python3

import os
import subprocess
import sys
import time


def verify(input_filename, output_filename, time = None):
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
                print("Error: Flight on output is not in input: {} {} {} {}".format(f, t, d, p))
                return 0
            calculated_price += input_flights[(f, t, d)]
            output_flights.append((f, t, d))
        output.close()
    except:
        print("Error: Wrong input -- error on the very first line")
        return 0

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

    print("OK -- price: {}, points: {:.2f}, weighted: {:.2f}{}".format(
        output_price,
        standard_price * 100.0 / output_price,
        standard_price * 100.0 / output_price * weight,
        "" if not time else ", time: {:.3f}".format(time),
    ))
    return standard_price * 100.0 / output_price * weight


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
        score += verify(fin, fout, time_end - time_start)
    print("\nSum: {:.2f}".format(score))
elif len(sys.argv) == 3:
    verify(sys.argv[1], sys.argv[2])
else:
    print("Expecting\neither 1 argument - Input directory,\nor 2 arguments - 1) Input file, 2) Contestant\'s output file")
    sys.exit(1)
