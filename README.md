### Kiwi challenge 2018

Main code - technical details:
- C++
- ugly
- parts adapted from last year's attempt in C
- random generator mt19937
  * tried also its several times faster variant SFMT, but dropped it as SIMD is not available
- data structures:
  * 3-dimensional array of prices, with days being the least significant index (to optimize caching)
  * array of "siblings" to a given airport, to quickly find airports in the same region

Algorithm:
- route = sequence of airports (not e.g. regions, or days)
- greedy initialization (a few variants of it, leave just the best)
- available reconfiguration moves (how to change the route):
  * swap two airports
  * swap three airports
  * Lin 2-opt
  * Lin 3-opt
  * replacement with another airport in the same region
  * cut an airport, paste it elsewhere in the route
- when applying the move, compute price difference rather than evaluate the new route
- get to local minimum: steepest descent using the moves
- simulated annealing using the moves
  * single temperature (tried one for each move but without much success)
  * fixed starting temperature 300
  * cooling schedule: T = 0.99 * T after X steps, or after X moves have been accepted by the annealing
  * simple bandit to prefer move types which have lead to a better price
  * "pressure": unavailable flights initially have price = 10 * maximal flight price, then geometric growth to 999999
- get to local minimum again: steepest descent using the moves

Support (Python):
- evaluation and verification, adapted from the official last year's verification script
- code for a simple search for the best parameters
- a few large semi-artificial input files
- a few ridiculously tiny input files
