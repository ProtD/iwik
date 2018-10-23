
#define DEBUG (1)
//#define NO_TIME_LIMIT (1)
//#define NDEBUG


#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <cstdarg>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
using namespace std;


// Sem vlozit sfmt.cpp, pokud je k dispozici SIMD (coz vypada, ze neni).

mt19937::result_type seed = time(0);

//auto random_generator = //sfmt19937(seed);
auto random_generator = mt19937(seed);

auto rand_real = std::bind(std::uniform_real_distribution<double>(0,1), random_generator);
auto rand_integer = std::bind(std::uniform_int_distribution<unsigned int>(0,0x8ffffffe), random_generator);
function<unsigned int()> rand_one_to_reg_minus_one_aux = std::bind(std::uniform_int_distribution<unsigned int>(0,1), random_generator);
function<unsigned int()> rand_one_to_reg_aux = std::bind(std::uniform_int_distribution<unsigned int>(0,1), random_generator);
inline int rand_int(int less_than) { return rand_integer() % less_than; }



#ifdef DEBUG
#include <signal.h>
#include <unistd.h>
#endif

#define MAXN (300+1)
#define MAX_PRICE (999999)
#define CODE_FIRST_CHAR ('0')
#define CODE_CHAR_SIZE ('z'-'0'+1)


int city_count = 0;
int region_count;
char city_codes[MAXN+1][3+1];
int initial;
int time_limit;
float wandering_started;
float wandering_interval;
int signal_caught = 0;
int pressure_price = MAX_PRICE;
int MAX_INT = std::numeric_limits<int>::max();

int two_way_city_count;
int city_stats[MAXN][2];   // ke kazdemu mestu, jestli ma vubec nejake vstupni / vystupni lety
#define INWARD (0)
#define OUTWARD (1)
int stats_min_price = 0, stats_max_price = 0, stats_mean_price = 0;  // let
int stats_price_count = 0, stats_price_max = 0;    // cela cesta
float stats_price_mean, stats_price_variance;
long stats_price_sum = 0;
long stats_price_square_dev_sum = 0;
struct {
  int min, max, sum, count;
} day_stats[MAXN];


int cached_1 = 0;
inline int rand_one_to_reg_minus_one() {
  if (cached_1 > 0) {
    int x = cached_1;
    cached_1 = 0;
    return x;
  }
  cached_1 = rand_one_to_reg_minus_one_aux();
  int x = 1 + cached_1 % (region_count-1);
  cached_1 = 1 + cached_1 / (region_count-1);
  return x;
}
int cached_2 = 0;
inline int rand_one_to_reg() {
  if (cached_2 > 0) {
    int x = cached_2;
    cached_2 = 0;
    return x;
  }
  cached_2 = rand_one_to_reg_aux();
  int x = 1 + cached_2 % region_count;
  cached_2 = 1 + cached_2 / region_count;
  return x;
}

class Siblings {
  int brothers[MAXN][MAXN]; // ke kazdemu mestu seznam mest ve stejnem regionu, prvni je pocet, konci 0
public:
  // vytvori z mest od i az do j sourozence
  inline int set_interval(int first_i, int last_i) {
    for (int j=first_i; j<=last_i; j++) {
      brothers[j][0] = last_i - first_i;
      for (int k=first_i; k<j; k++)
        brothers[j][k - first_i + 1] = k;
      for (int k=j+1; k<=last_i; k++)
        brothers[j][k - first_i] = k;
      brothers[j][last_i + 1 - first_i] = 0;
    }
  }

  // vrati order-teho sourozence mesta i (nekde mezi nimi je i mesto samotne)
  inline int get(int i, int order) {
    assert(order < count(i));
    assert(order == 0 || brothers[i][order] > 0);
    return (order == 0) ? i : brothers[i][order];
  }

  // pocet sourozencu mesta i, vcetne jeho samotneho
  inline int count(int i) { return brothers[i][0] + 1; }

  // odstrani mesto i
  inline int remove(int i) {
    for (int k=1; k<count(i); k++) {
      int sibl = get(i, k);
      int found = 0;
      for (int j=1; j<count(i); j++) {
        if (get(sibl, j) == i)
          found++;
        else
          brothers[sibl][j-found] = brothers[sibl][j];
      }
      brothers[sibl][0] -= found;
    }
    brothers[i][0] = 0;
  }

  // vrati nahodneho sourozence mesta i
  inline int get_random_sibling(int i) {
    int cnt = count(i);
    if (cnt == 1)
      return 0;
    return brothers[i][1 + rand_int(cnt - 1)];
  }

  // vrati nahodneho sourozence mesta i, vcetne jeho sameho
  inline int get_random_sibling_or_me(int i) {
    int cnt = count(i);
    if (cnt == 1)
      return i;
    int rnd = rand_int(cnt);
    if (rnd == 0)
      return i;
    else
      return brothers[i][rnd];
  }
};
Siblings siblings;

class Prices {
public:
  int *matrix;  // 1. den je indexovan od 0, kdezto 2. mesto a 3. mesto jsou indexovany od 1
  int region_count_times_city_count_plus_one;
  
  Prices() {
    // // alokace m a inicializace maximalni cenou
    // for (int k=1; k<=city_count; k++)
    //    m[0][1][k] = MAX_PRICE;
    // for (int i=0; i<city_count; i++) // usetrime par milisekund, kdyz to nebudem inicializovat cely
    //  for (int j=1; j<=city_count; j++)
    //    std::copy(m[0][1], m[0][1] + city_count, m[i][j]);

    region_count_times_city_count_plus_one = region_count * (city_count+1);
    matrix = new int[region_count * (city_count+1) * (city_count+1)];
    //todo: zrychlit?
    for (int from=1; from<=city_count; from++)
      for (int to=1; to<=city_count; to++)
        for (int day=0; day<region_count; day++)
          matrix[index(from, to, day)] = MAX_PRICE;
  }
  ~Prices() {
    delete[] matrix;
  }
  
  inline int index(int from, int to, int day) {
    assert((from)>0 && (from)<=city_count);
    assert((to)>0 && (to)<=city_count);
    assert((day)>=0 && (day)<region_count);
    return day + from * region_count + to * region_count_times_city_count_plus_one;
  }

  inline void set(int from, int to, int day, int price) {
    int p = matrix[index(from, to, day)];
    if (price > stats_price_max)
      stats_price_max = price;
    if (p == MAX_PRICE) {
      matrix[index(from, to, day)] = price;
      stats_price_sum += price;
      stats_price_count += 1;
    } else if (price < p) {
      matrix[index(from, to, day)] = price;
      stats_price_sum += price - p;
    }
  }
};

Prices *prices;
inline int PRICE(int from, int to, int day) {
  int x = prices->matrix[prices->index(from, to, day)];
  if (x == MAX_PRICE)
    return pressure_price;
  else
    return x;
}


// cas od pocatku v sekundach
struct timeval tv_start;
float elapsed()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return (tv.tv_sec - tv_start.tv_sec) + (tv.tv_usec - tv_start.tv_usec) / 1000000.0f;
}

// kdyz uz je mesto zname, vrati jeho ciselne oznaceni
// kdyz jeste neni zname, prida ho mezi zname, a vrati 0
int a2i(int *city_list, char *code) {
  int i = (code[0] - CODE_FIRST_CHAR) * CODE_CHAR_SIZE * CODE_CHAR_SIZE + (code[1] - CODE_FIRST_CHAR) * CODE_CHAR_SIZE + (code[2] - CODE_FIRST_CHAR);
  if (city_list[i] > 0)
    return city_list[i];
  city_list[i] = ++city_count;
  strncpy(city_codes[city_count], code, 3);
  city_codes[city_count][4] = '\0';
  return city_count;
}

void load() {
  int *city_list = (int *) calloc(CODE_CHAR_SIZE * CODE_CHAR_SIZE * CODE_CHAR_SIZE, sizeof(int));
  const size_t line_size = 300;
  char *line = (char *) malloc(line_size);
  char *c;
  char initial_code[3+1] = "";

  // uvodni radek
  fgets(line, line_size, stdin);
  region_count = 0;
  for (c=line; c[0] != ' '; c++)
    region_count = region_count*10 + (c[0]-'0');
  strncpy(initial_code, ++c, 3);

  // definice oblasti
  int region_from = 1;
  for (int i=0; i<region_count; i++) {
    fgets(line, line_size, stdin);
    fgets(line+1, line_size, stdin);
    c = line;
    int region_to = region_from;
    do {
      c++;
      a2i(city_list, c);
      c += 3;
      region_to++;
    } while (c[0] != '\n');
    siblings.set_interval(region_from, region_to-1);
    region_from = region_to;
  }

  if (region_count <= 20 && city_count < 50)
    time_limit = 3;
  else if (region_count <= 100 && city_count < 200)
    time_limit = 5;
  else
    time_limit = 15;

#ifdef DEBUG
  printf("%#6.3fs Nacteny oblasti.\n", elapsed());
#endif

  prices = new Prices();
  memset(city_stats, 0, MAXN*2*sizeof(int));
  
#ifdef DEBUG
  printf("%#6.3fs Inicializovana matice cen.\n", elapsed());
#endif

  // nacteni cen
  while (fgets(line, line_size, stdin) != NULL) {
    int from = a2i(city_list, line);
    int to = a2i(city_list, line+4);

    int day = 0;
    for (c = line+8; c[0] != ' '; c++)
      day = day*10 + (c[0]-'0');

    int price = 0;
    for (c++; c[0] != '\n' && c[0] != '\0'; c++)
      price = price*10 + (c[0]-'0');

    if (day > 0)
      prices->set(from, to, day-1, price);
    else
      for (int j=0; j<region_count; j++)
        prices->set(from, to, j, price);
    city_stats[from][OUTWARD] = 1;
    city_stats[to][INWARD] = 1;
  }
  free(line);

  initial = a2i(city_list, initial_code);
  for (int order=0; order<siblings.count(initial); order++) {
    int k = siblings.get(initial, order);
    if (k == initial)
      city_stats[k][INWARD] = 1;
    else 
      city_stats[k][OUTWARD] = 1;
  }

#ifdef DEBUG
  cout << "Nalezene slepe ulicky:";
#endif
  two_way_city_count = city_count;
  for (int i=1; i<=city_count; i++)
    if (city_stats[i][INWARD] == 0 || city_stats[i][OUTWARD] == 0) {
      siblings.remove(i);
      two_way_city_count--;
#ifdef DEBUG
      cout << " " << city_codes[i];
    }
  cout << " ... zbyva " << two_way_city_count << endl;
#else
}
#endif

stats_price_mean = float(stats_price_sum) / stats_price_count;
for (int day=0; day<region_count; day++) {
  int *min = &day_stats[day].min;
  int *max = &day_stats[day].max;
  int *sum = &day_stats[day].sum;
  int *count = &day_stats[day].count;
  *min = MAX_PRICE;
  *max = *sum = *count = 0;
  for (int j=1; j<=city_count; j++)
    for (int k=1; k<=city_count; k++) {
      int p = PRICE(j, k, day);
      if (p < MAX_PRICE) {
        stats_price_square_dev_sum += (stats_price_mean - p) * (stats_price_mean - p);
        *sum += p;
        if (p > *max)
          *max = p;
        if (p < *min)
          *min = p;
        (*count)++;
      }
    }
 }
stats_price_variance = sqrt(float(stats_price_square_dev_sum) / stats_price_count);

for (int i=0; i<region_count; i++) {
  stats_min_price += day_stats[i].min;
  stats_max_price += day_stats[i].max;
  stats_mean_price += day_stats[i].sum / day_stats[i].count;
 }
#ifdef DEBUG
cout << "Statistiky:" << endl
<< "Cena minimalni / maximalni / prumerne cesty: " << stats_min_price << ", " << stats_max_price << ", " << stats_mean_price << endl
<< "Pocet spoju: " << stats_price_count << endl
<< "Prumerna / maximalni cena spoje: " << stats_price_mean << ", " << stats_price_max << endl
<< "Rozptyl: " << stats_price_variance << endl;
  
printf("%#6.3fs Nacteno %d mest, %d oblasti, limit je %d s, zacina se v %s (%d).\n", elapsed(), city_count, region_count, time_limit, initial_code, initial);
#endif

free(city_list);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class Route {
public:
  int *route;
  
  Route() { route = new int[region_count+1]; }
  ~Route() { delete route; }
  
  inline void copy(const Route &r2) {
    std::copy(r2.route, r2.route + region_count + 1, route);
  }

  inline Route(const Route& r2) {
    route = new int[region_count+1];
    Route::copy(r2);
  }
  
  inline int price(int from = 1, int to = region_count) {
    int price = 0;
    for (int i=from-1; i<to; i++)
      price += PRICE(route[i], route[i+1], i);
    return price;
  }

  inline void print();

  void random_init();
  void greedy_init(int forward, int backward, Route *anti);

  // prohozeni i a j
  inline void two_swap(int i, int j);
  inline int two_points_price(int i, int j, int a, int b);
  inline int price_after_two_swap(int i, int j, int current_price);
  int best_two_swap(int p = 0);
  int gd_two_swap(int p = 0, int until = 0);

  // rotace i, j a k
  inline void three_swap(int i1, int i2, int i3);
  inline int three_points_price(int i, int j, int k, int a, int b, int c);
  inline int price_after_three_swap(int i, int j, int k, int current_price);

  // nahrazeni i mestem j (ze stejneho regionu)
  inline void brother_opt(int i, int k);
  inline int price_after_brother_opt(int i, int brother, int current_price);
  int best_brother(int p = 0);

  // rotace useku mezi i a j, predpoklada i < j
  inline void two_opt(int i1, int j2);
  inline int price_after_two_opt(int i1, int i2, int current_price, int max);
  inline int compute_two_opt(int i1, int i2, int max, int only_changes);
  int best_two_opt(int p = 0);
  int gd_two_opt(int p = 0);

  // prohozeni/rotace tri useku, predpoklada i < j < k
  inline void three_opt(int which, int i1, int i2, int i3);
  inline int price_after_three_opt(int which, int i1, int i2, int i3, int current_price, int max);
  int best_three_opt(int proximity, int p = 0);
  int gd_three_opt(int proximity, int p = 0);
};
  
Route *result_route = NULL;
int result_price = MAX_INT;

inline int soft_limit_passed() {
  if (result_price < pressure_price * region_count)
    return elapsed() + region_count / 500.0 + 0.1 > time_limit;
  else
    return elapsed() + 0.1 > time_limit;
}

inline int hard_limit_passed() {
  return elapsed() + 0.05 > time_limit;
}

void Route::print() {
  cout << price() << endl;
  for (int i=0; i<region_count; i++)
    cout << city_codes[route[i]] << " "
         << city_codes[route[i+1]] << " "
         << i+1 << " "
         << PRICE(route[i], route[i+1], i)
         << endl << flush;
}

// nahodna inicializace
void Route::random_init() {
  // inicializace poporade
  route[0] = initial;
  int j = 1;
  for (int i=0; i<region_count; i++) {
    if (j < initial)
      route[i+1] = siblings.get_random_sibling_or_me(j);
    else if (j >= initial + siblings.count(initial))
      route[i] = siblings.get_random_sibling_or_me(j);
    j += siblings.count(j);
  }
  route[region_count] = siblings.get_random_sibling_or_me(initial);

  // nahodne zprehazeni
  for (int i=1; i<region_count; i++) {
    int j = rand_int(region_count-i);
    swap(route[region_count-1-j], route[i]);
  }
}

// greedy inicializace -- pridava na jeden z koncu vzdy nejlevnejsi usek
// forward = pridava na konec, backward = pridava na zacatek, muze byt zapnute oboje
// anti = cesta, ktere by to melo byt pokud mozno nepodobne
void Route::greedy_init(int forward, int backward, Route *anti) {
  route[0] = initial;
  route[region_count] = initial;

  int was[MAXN+1];
  memset(was, 0, MAXN*sizeof(int));
  for (int order=0; order<siblings.count(initial); order++)
    was[siblings.get(initial, order)] = 1;
  for (int i=1; i<=city_count; i++)
    if (city_stats[i][INWARD] == 0 || city_stats[i][OUTWARD] == 0)
      was[i] = 1;

  int last[2] = { initial, initial };
  int is[2] = { 1, region_count-1 };
  while (1) {
    int min_price = MAX_PRICE;
    int min_j = 0;
    int min_direction = -1;
    for (int direction = 0; direction <= 1; direction++) {
      if (direction == 0 && !forward)
        continue;
      if (direction == 1 && !backward)
        continue;
      int discouraged_a = -1;
      int discouraged_b = -1;
      if (anti != NULL)
        for (int j=0; j<region_count; j++)
          if (anti->route[j] == last[direction]) {
            discouraged_a = anti->route[(j-1) % region_count];
            discouraged_b = anti->route[(j+1) % region_count];
          }
      for (int j=1; j<=city_count; j++)
        if (!was[j]) {
          min_j = j;
          min_direction = direction;
        }
      for (int j=1; j<=city_count; j++)
        if (!was[j]) {
          int price = (direction == 0) ? PRICE(last[0], j, is[0]-1) : PRICE(j, last[1], is[1]);
          if (j == discouraged_a || j == discouraged_b)
            price += MAX_PRICE / 10;
          if (price < min_price) {
            /*
              #ifdef DEBUG
              if (direction == 0)
              printf("p %d -> %d (%d): %d\n", last[0], j, is[0]-1, price);
              else
              printf("p %d -> %d (%d): %d\n", j, last[1], is[1], price);
              #endif
            */
            min_price = price;
            min_j = j;
            min_direction = direction;
          }
        }
    }
    route[is[min_direction]] = min_j;
    for (int order=0; order<siblings.count(min_j); order++)
      was[siblings.get(min_j, order)] = 1;
    last[min_direction] = min_j;
    if (min_direction == 0) {
      is[0] += 1;
    } else {
      is[1] -= 1;
    }
    if (is[0] > is[1])
      break;
  }
}

void Route::two_swap(int i, int j) {
  swap(route[i], route[j]);
}

int Route::two_points_price(int i, int j, int a, int b) {
  if (j == i + 1)
    return 
      + PRICE(route[i-1], route[a], i-1)
      + PRICE(route[a], route[b], i)
      + PRICE(route[b], route[j+1], j);
  else if (j == i - 1)
    return 
      + PRICE(route[j-1], route[b], j-1)
      + PRICE(route[b], route[a], j)
      + PRICE(route[a], route[i+1], i);
  else
    return
      + PRICE(route[i-1], route[a], i-1)
      + PRICE(route[a], route[i+1], i)
      + PRICE(route[j-1], route[b], j-1)
      + PRICE(route[b], route[j+1], j);
}

int Route::price_after_two_swap(int i, int j, int current_price) {
  return current_price - two_points_price(i, j, i, j) + two_points_price(i, j, j, i);
}

int Route::best_two_swap(int p) {
  Route minr = Route();
  if (p == 0)
    p = price();
  int minp = p;
  for (int i=1; i<region_count; i++)
    for (int j=i+1; j<region_count; j++) {
      int p2 = price_after_two_swap(i, j, p);
      if (p2 < minp) {
        minp = p2;
        minr.copy(*this);
        minr.two_swap(i, j);
        assert(minr.price() == minp);
      }
    }
  if (minp < p)
    Route::copy(minr);
  return minp;
}

int Route::gd_two_swap(int p, int until) {
  if (p == 0)
    p = price();
  int prevp;
  do {
    prevp = p;
    p = best_two_swap(p);
  } while (p < prevp && p > until && !hard_limit_passed());
  return p;
}

void Route::three_swap(int i, int j, int k) {
  int x = route[i];
  route[i] = route[j];
  route[j] = route[k];
  route[k] = x;
}

// pokus o spocitani ceny po three-swap - ale bohuzel zcela nefunguje
int Route::three_points_price(int i, int j, int k, int a, int b, int c) {
  return
    + PRICE(route[i-1], route[a], i-1)
    + PRICE(route[a], route[i+1], i)
    + PRICE(route[j-1], route[b], j-1)
    + PRICE(route[b], route[j+1], j)
    + PRICE(route[k-1], route[c], k-1)
    + PRICE(route[c], route[k+1], k);
}

int Route::price_after_three_swap(int i, int j, int k, int current_price) {
  if (((i-j>1) || (i-j<-1)) && ((j-k>1) || (j-k<-1)) && ((k-i>1) || (k-i<-1))) {
    return current_price - three_points_price(i, j, k, i, j, k) + three_points_price(i, j, k, j, k, i);
  } else {
    // moc slozity, tak to holt napocitame cele :-)
    Route r2 = Route(*this);
    r2.three_swap(i, j, k);
    return r2.price();
  }
}
  
void Route::brother_opt(int i, int k) {
  route[i] = k;
}

int Route::price_after_brother_opt(int i, int k, int current_price) {
  int result = current_price
    - PRICE(route[i-1], route[i], i-1)
    + PRICE(route[i-1], k, i-1);
  if (i < region_count)
    return result
      - PRICE(route[i], route[i+1], i)
      + PRICE(k, route[i+1], i);
  return result;  
}

int Route::best_brother(int p) {
  Route minr = Route();
  if (p == 0)
    p = price();
  int minp = p;
  for (int i=1; i<region_count; i++)
    for (int j=0; j<siblings.count(route[i]); j++) {
      int k = siblings.get(route[i], j);
      if (k == i)
        continue;
      int p2 = price_after_brother_opt(i, k, p);
      if (p2 < minp) {
        minr.copy(*this);
        minr.brother_opt(i, k);
        minp = p2;
      }
    }
  if (minp < p)
    Route::copy(minr);
  return minp;
}

void Route::two_opt(int i1, int i2) {
  for (int k=i1, l=i2-1; k < l; k++, l--)
    swap(route[k], route[l]);
}

int Route::price_after_two_opt(int i1, int i2, int current_price, int max) {
  if (i2-i1 < region_count/2 && current_price > 0) {
    int changes = compute_two_opt(i1, i2, MAX_INT, 1);
    if (changes == MAX_INT)
      return max;
    else
      return current_price - price(i1, i2) + changes;
  }
  else
    return compute_two_opt(i1, i2, max, 0);
}

int Route::compute_two_opt(int i1, int i2, int max, int only_changes) {
  int j, day;
  int price = 0;
  if (!only_changes) {
    for (j=0; j<i1-1; j++) {
      price += PRICE(route[j], route[j+1], j);
    }
    if (price > max)
      return max;
  }
  price += PRICE(route[i1-1], route[i2-1], i1-1);
  for (j=i2-1, day=i1; j>i1; j--, day++) {
    price += PRICE(route[j], route[j-1], day);
  }
  price += PRICE(route[i1], route[i2], i2-1);
  if (!only_changes) {
    if (price > max)
      return max;
    for (j=i2; j<region_count; j++)
      price += PRICE(route[j], route[j+1], j);
  }
  return price;
}

int Route::best_two_opt(int p) {
  Route minr = Route();
  if (p == 0)
    p = price();
  int minp = p;
  for (int i=1; i<region_count; i++)
    for (int j=i+1; j<=region_count; j++) {
      int p2 = price_after_two_opt(i, j, p, minp);
      if (p2 < minp) {
        minp = p2;
        minr.copy(*this);
        minr.two_opt(i, j);
        assert(minr.price() == minp);
      }
    }
  if (minp < p)
    Route::copy(minr);
  return minp;
}

int Route::gd_two_opt(int p) {
  if (p == 0)
    p = price();
  int prevp;
  do {
    prevp = p;
    p = best_two_opt(p);
  } while (p < prevp && !hard_limit_passed());
  return p;
}

void Route::three_opt(int which, int i1, int i2, int i3) {
  switch (which) {
  case 0: {
    for (int k=i1, l=i2-1; k < l; k++, l--)
      swap(route[k], route[l]);
    for (int k=i2, l=i3-1; k < l; k++, l--)
      swap(route[k], route[l]);
    break;
  }
  case 1: {
    if (i3-i2 > i2-i1) {
      int foo[i3-i2];
      memcpy(foo, route + i2, (i3-i2) * sizeof(int));
      memcpy(route + i3 - (i2-i1), route + i1, (i2-i1) * sizeof(int));
      memcpy(route + i1, foo, (i3-i2) * sizeof(int));
    } else {
      int foo[i2-i1];
      memcpy(foo, route + i1, (i2-i1) * sizeof(int));
      memcpy(route + i1, route + i2, (i3-i2) * sizeof(int));
      memcpy(route + i3 - (i2-i1), foo, (i2-i1) * sizeof(int));
    }
    break;
  }
  }
}

// nefunguje??
int Route::price_after_three_opt(int which, int i1, int i2, int i3, int current_price, int max) {
  int j, day;
  int p = 0;
  if (current_price > 0)
    p = current_price - price(i1, i3);
  switch (which) {
  case 0: {
    p += PRICE(route[i1-1], route[i2-1], i1-1);
    for (j=i2-1, day=i1; j>i1; j--, day++) {
      p += PRICE(route[j], route[j-1], day);
      if (p > max)
        return max;
    }
    p += PRICE(route[i1], route[i3-1], i2-1);
    for (j=i3-1, day=i2; j>i2; j--, day++) {
      p += PRICE(route[j], route[j-1], day);
      if (p > max)
        return max;
    }
    p += PRICE(route[i2], route[i3], i3-1);
    break;
  }
  case 1: {
    p += PRICE(route[i1-1], route[i2], i1-1);
    for (j=i2, day=i1; j<i3-1; j++, day++) {
      p += PRICE(route[j], route[j+1], day);
      if (p > max)
        return max;
    }
    p += PRICE(route[i3-1], route[i1], i1+i3-i2-1);
    for (j=i1, day=i1+i3-i2; j<i2-1; j++, day++) {
      p += PRICE(route[j], route[j+1], day);
      if (p > max)
        return max;
    }
    p += PRICE(route[i2-1], route[i3], i3-1);
    break;
  }
  }
  return p;
}
                                        
int Route::best_three_opt(int proximity, int p) {
  Route minr = Route();
  if (p == 0)
    p = price();
  int minp = p;
  for (int i=1; i<region_count-1; i++) {
    int bound = i + proximity;
    if (bound > region_count)
      bound = region_count;
    for (int k=i+2; k<bound; k++) {
      int part_price = price(i, k);
      int limit_max = minp + part_price - p;
      for (int j=i+1; j<k; j++)
        for (int which=0; which<=1; which++) {
          int p2 = p - part_price + price_after_three_opt(which, i, j, k, 0, limit_max);
          //int p2 = price_after_three_opt(which, i, j, k, p, minp);       // tohle by melo spocitat to same, ale pomaleji
          if (p2 < minp) {
            minp = p2;
            limit_max = minp + part_price - p;
            minr.copy(*this);
            minr.three_opt(which, i, j, k);
            assert(minr.price() == minp);
          }
        }
    }
  }
  if (minp < p)
    Route::copy(minr);
  return minp;
}

int Route::gd_three_opt(int proximity, int p) {
  if (p == 0)
    p = price();
  int prevp;
  do {
    prevp = p;
    p = best_three_opt(proximity, p);
  } while (p < prevp && !hard_limit_passed());
  return p;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void checkout(Route &route, int price) {
  if (price < result_price) {
    result_route->copy(route);
    result_price = price;
#ifdef DEBUG
    printf("%#6.3fs ----------------------- Nalezena cena: %d.\n", elapsed(), result_price);
#endif
  }
}

#ifdef DEBUG
void signal_handler(int i){
  printf("Zachycen signal %d.\n", i);
  signal_caught++;
  if (signal_caught > 5)
    exit(1);
}
#endif

const int st_none = 0;
const int st_two_swap = 1;
const int st_three_swap = 2;
const int st_two_opt = 3;
const int st_three_opt = 4;
const int st_brother_opt = 5;
const int st_size = 6;

int random_choice(int count, int sum, ...) {
  va_list args;
  va_start(args, sum);
  int dice = rand_int(sum);
  int so_far = 0;
  for (int i=0; i<count; i++) {
    so_far += va_arg(args, int);
    int step_type = va_arg(args, int);
    if (dice < so_far)
      return step_type;
  }
  va_end(args);
  return st_none;
}

enum step_result_t { sr_forbidden, sr_rejected, sr_accepted, sr_descended, sr_enhanced };

inline step_result_t wander_step(Route &r, int &p, Route &bestr, int &bestp,
                                 const float &normalized_temperature, const unsigned int step,
                                 const int *step_type_probability, int &step_type) {
  int p2;
  int step_params[4];

  step_type = st_none;
  do {
    step_type = random_choice(5, step_type_probability[0],
                              step_type_probability[st_two_swap],    st_two_swap,
                              step_type_probability[st_three_swap],  st_three_swap,
                              step_type_probability[st_two_opt],     st_two_opt,
                              step_type_probability[st_three_opt],   st_three_opt,
                              step_type_probability[st_brother_opt], st_brother_opt);
    switch (step_type) {

    case st_brother_opt:
      int to_change;
      do {
        step_params[0] = rand_one_to_reg();
        to_change = siblings.get_random_sibling(r.route[step_params[0]]);
      } while (to_change == 0);
      step_params[1] = to_change;
      p2 = r.price_after_brother_opt(step_params[0], to_change, p);
      break;
      
    case st_two_swap:
    case st_three_swap:
    case st_two_opt:
    case st_three_opt:
      int indices = (step_type == st_three_swap || step_type == st_three_opt) ? 3 : 2;
      for (int i=0; i<indices; i++) {
        int leave;
        do {
          step_params[i] = rand_one_to_reg_minus_one();
          leave = 1;
          for (int j=0; j<i; j++)
            if (step_params[i] == step_params[j])
              leave = 0;
        } while (!leave);
      }
      if (step_type != st_two_swap && step_type != st_three_swap) {
        sort(step_params, step_params+indices);
        if (step_params[indices-1] - step_params[0] > 1000) {       // todo
          step_type = st_none;
          break;
        }
      }
        
      switch (step_type) {
      case st_two_swap:   p2 = r.price_after_two_swap(step_params[0], step_params[1], p); break;
      case st_three_swap: p2 = r.price_after_three_swap(step_params[0], step_params[1], step_params[2], p); break;
      case st_two_opt:    p2 = r.price_after_two_opt(step_params[0], step_params[1], p, p + MAX_PRICE); break;
      case st_three_opt:
        step_params[3] = rand_int(2);
        p2 = r.price_after_three_opt(step_params[3], step_params[0], step_params[1], step_params[2], p, p + MAX_PRICE);
        break;
      }
    }
  } while (step_type == st_none);

  if (p2 - p >= MAX_PRICE)
    return sr_forbidden;

  int ascended = (p2 > p);
  int accepted;
  if (ascended)
    accepted = (rand_real() < exp((p-p2)/normalized_temperature));
  else
    accepted = 1;

  if (p2 < bestp || accepted) {
    Route *r2;
    if (accepted) {
      p = p2;
      r2 = &r;  // nemusime nikam kopirovat
    } else {
      r2 = new Route(r);
    }
    
    switch (step_type) {
    case st_brother_opt: r2->brother_opt(step_params[0], step_params[1]); break;
    case st_two_swap:    r2->two_swap(step_params[0], step_params[1]); break;
    case st_three_swap:  r2->three_swap(step_params[0], step_params[1], step_params[2]); break;
    case st_two_opt:     r2->two_opt(step_params[0], step_params[1]); break;
    case st_three_opt:   r2->three_opt(step_params[3], step_params[0], step_params[1], step_params[2]); break;
    }
    if (r2->price() != p2)
      cout << step << ": typ " << step_type << ": stated price " << p2 << " != real " << r2->price() << " (diff. " << p2-r2->price() << "), "
           << step_params[0] << ", " << step_params[1] << ", " << step_params[2] << endl;
    assert(r2->price() == p2);
      
    if (p2 < bestp) {
      bestp = p2;
      bestr.copy(*r2);
      checkout(*r2, p2);
      if (!accepted)
        delete r2;
      return sr_enhanced;
    }
    if (ascended)
      return sr_accepted;
    else if (accepted)
      return sr_descended;
  }
  return sr_rejected;
}

inline float exp_interpolate(float fraction, float from, float to) {
  return from * exp(fraction * log(to / from));
}

inline float update_pressure_price_get_temperature(float starting_temperature) {
  float current_time = elapsed();
  if (current_time > time_limit)
    pressure_price = MAX_PRICE;
  else
    pressure_price = exp_interpolate((current_time - wandering_started) / wandering_interval,
                                     2 * stats_price_max, MAX_PRICE);
  if (pressure_price > MAX_PRICE)
    pressure_price = MAX_PRICE;
  result_price = result_route->price();

  current_time -= floor(current_time / time_limit) * time_limit;
  return exp_interpolate((current_time - wandering_started) / wandering_interval,
                         //                       stats_price_mean / 5, stats_price_mean / 10);
                         region_count / 2, region_count / 8);
}

int wander(Route *rs, int route_count, float starting_temperature) {
  int ps[route_count];
  int bestps[route_count];
  Route bestrs[route_count];
  int step_types_cnt[st_size] = {0, 1, 1, 1, 1, 1};
  int step_types_enhanced[st_size] = {0, 10, 4, 2, 1, 3};
  int step_types_accepted[st_size] = {0, 0, 0, 0, 0, 0};
  int step_types_probability[st_size];

  wandering_started = elapsed();
  wandering_interval = (time_limit - (wandering_started + region_count / 500.0 + 0.1));
  
  float temperature = update_pressure_price_get_temperature(starting_temperature);
  for (int i=0; i<route_count; i++) {
    bestrs[i].copy(rs[i]);
    bestps[i] = ps[i] = rs[i].price();
    checkout(rs[i], ps[i]);
  }

  unsigned int step = 0;
  do {
#ifdef DEBUG
    if (step % 1000000 == 0) {
      printf("%#6.3fs Krok %d M, teplota %.4f, tlak %d, akt.cena %d, nej.cena %d, bandita: ",
             elapsed(), step / 1000000, temperature, pressure_price, ps[0], result_price);
      for (int i=1; i<st_size; i++)
        cout << i << ":" << step_types_enhanced[i] << "/" << step_types_accepted[i] << "/" << step_types_cnt[i] / 1000 << "k ";
      cout << endl;
    }
#endif
    if (step % 10000 == 0) {
      step_types_probability[0] = 0;
      for (int i=1; i<st_size; i++) {
        if (i == st_brother_opt && two_way_city_count == region_count)
          step_types_probability[i] = 0;
        else
          step_types_probability[i] = ceil(step_types_enhanced[i] * 100000.0 / step_types_cnt[i]);
        step_types_probability[0] += step_types_probability[i];
      }
    }

    for (int i=0; i<route_count; i++) {
      int step_type;
      step_result_t step_result = wander_step(rs[i], ps[i], bestrs[i], bestps[i], temperature, step, step_types_probability, step_type);
      step_types_cnt[step_type]++;
      if (step_result == sr_enhanced)
        step_types_enhanced[step_type]++;
      if (step_result == sr_accepted || step_result == sr_descended || step_result == sr_enhanced)
        step_types_accepted[step_type]++;
    }

    step++;
    if (step % 10000 == 0) {
      temperature = update_pressure_price_get_temperature(starting_temperature);
      for (int i=0; i<route_count; i++) {
        ps[i] = rs[i].price();
        bestps[i] = bestrs[i].price();
      }
    }

    // trocha genetiky: naklonujeme lepsi polovinu na ukor horsi poloviny
    if (step % 100000 == -1) {
    outer:
      for (int i=0; i<route_count; i++) {
        for (int j=0; j<i; j++)
          if (bestps[i] < bestps[j]) {
            swap(rs[i].route, rs[j].route);
            swap(ps[i], ps[j]);
            swap(bestps[i], bestps[j]);
            swap(bestrs[i], bestrs[j]);
            goto outer;
          }
      }
      int last_bestp = 0;
      int j = route_count - 1;
      for (int i=0; i<j; i++) {
        if (bestps[i] == last_bestp)  // neklonujeme ty se stejnou cenou -- pokus o diverzitu
          continue;
        rs[j].copy(rs[i]);
        ps[j] = ps[i];
        bestps[j] = bestps[i];
        bestrs[j].copy(bestrs[i]);
        j--;
      }
    }

#ifdef NO_TIME_LIMIT
  } while (!signal_caught);
#elif defined DEBUG
} while ((step % 100 != 0 || !soft_limit_passed()) && !signal_caught);
#else
} while (step % 100 != 0 || !soft_limit_passed());
#endif
}

// moc nefunguje - o nekolik radu teplotu nadsazuje
float probe_temperature(Route *rs, int route_count)
{
  int p, bestp;
  Route r, bestr;
  int step_types_probability[st_size] = {1, 1, 0, 0, 0, 0};

  int orig_p[route_count];
  for (int i=0; i<route_count; i++)
    orig_p[i] = rs[i].price();

  float temperature = 100;
  int accepts, rejects;
  do {
    accepts = 0;
    rejects = 0;
    for (unsigned int step=0; step<10000; step++) {
      for (int i=0; i<route_count; i++) {
        int step_type;
        r.copy(rs[i]);
        p = orig_p[i];
        bestp = orig_p[i];
        step_result_t step_result = wander_step(r, p, bestr, bestp, temperature, step, step_types_probability, step_type);
        if (step_result == sr_accepted)
          accepts += 1;
        else if (step_result == sr_rejected)
          rejects += 1;
      }
    }
    cout << "Teplota " << temperature << ": " << accepts << " : " << rejects << endl;
    temperature *= 1.2;
  } while (accepts < rejects);
  return temperature;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {
  gettimeofday(&tv_start, 0);
  
#ifdef DEBUG
  struct sigaction sigIntHandler;
  sigIntHandler.sa_handler = signal_handler;
  sigemptyset(&sigIntHandler.sa_mask);
  sigIntHandler.sa_flags = 0;
  sigaction(SIGINT, &sigIntHandler, NULL);
#endif

  // inicializace promennych a nacteni vstupu
  load();

  // kdyz nechci zvedat skore
  //if (region_count == 10)
  //  exit(0);

  rand_one_to_reg_minus_one_aux = std::bind(std::uniform_int_distribution<unsigned int>(0,(region_count-1)*(region_count-1)-1), random_generator);
  rand_one_to_reg_aux = std::bind(std::uniform_int_distribution<unsigned int>(0,region_count*region_count-1), random_generator);

  result_route = new Route();

  int route_count = 6;
  Route rs[route_count];
  int ps[route_count];

  // zjisteni spravne pocatecni teploty
  //rs[0].random_init();
  //float starting_temperature = probe_temperature(rs, 1);
  //cout << "Spravna teplota: " << starting_temperature << endl;
  
  // inicializace cest
  rs[0].greedy_init(1, 0, NULL);
  rs[1].greedy_init(0, 1, NULL);
  rs[2].greedy_init(1, 1, NULL);
  rs[3].greedy_init(1, 0, &rs[0]);
  rs[4].greedy_init(0, 1, &rs[1]);
  rs[5].greedy_init(1, 1, &rs[2]);
  // for (int i=6; i<route_count; i++)
  //  rs[i].random_init();

  // sjedeme pokud mozno pod MAX_PRICE pomoci dva-swapu
  for (int i=0; i<route_count; i++) {
    ps[i] = rs[i].price();
#ifdef DEBUG
    printf("%#6.3fs %d inicializovano, cena: %d.\n", elapsed(), i, ps[i]);
#endif
    ps[i] = rs[i].gd_two_swap(ps[i], MAX_PRICE);
#ifdef DEBUG
    printf("%#6.3fs %d 2-swap, cena: %d.\n", elapsed(), i, ps[i]);
#endif
  }

  // nechame jen nejlepsi
  int best = 0;
  for (int i=0; i<route_count; i++) {
    ps[i] = rs[i].price();
    if (ps[i] < ps[best])
      best = i;
  }
  rs[0].copy(rs[best]);
  ps[0] = ps[best];
  route_count = 1;

  // sjedeme do minima pomoci dva- nebo tri-optu
  for (int i=0; i<route_count; i++) {
    int p = rs[i].best_brother();
    if (region_count < 100)
      p = rs[i].gd_three_opt(50);
    else
      p = rs[i].gd_two_opt();
#ifdef DEBUG
    printf("%#6.3fs %d 2/3-opt, cena: %d.\n", elapsed(), i, p);
#endif
    checkout(rs[i], p);
    if (soft_limit_passed()) {
#ifdef DEBUG
      printf("%#6.3fs dosel cas, koncime.\n", elapsed());
#endif
      break;
    }
  }

  // zihani
#ifdef DEBUG
  printf("%#6.3fs jdeme bloudit.\n", elapsed());
#endif
  wander(rs, route_count, 0);
  //wander(rs, route_count, 50);

  // vypis vysledku
#ifdef DEBUG
  printf("%#6.3fs Jeste dooptimalizujeme, zatim cena: %d.\n", elapsed(), result_price);
#endif
  if (region_count < 100)
    result_route->gd_three_opt(50);
  else
    result_route->gd_two_opt();
  result_route->print();

#ifdef DEBUG
  printf("%#6.3fs Konec, cena: %d.\n", elapsed(), result_route->price());
#endif

  delete prices;

  return 0;
}

/*
  Nejlepsi dosazene reseni / reseni za 100 bodu / vhodna teplota cca
  1:   1396 /   1396
  2:   1407 /   2159 /  20
  3:  39290 /  44151 /  60
  4: 115690 / 118814 / 120

  Statistiky:
  pr.  mest  oblasti     spoju   pr.cena  rozptyl
  1      10       10       900      1645     1921
  2     192       96    621120       130       45
  3     200      150   1702027       715      385
  4     300      300   3567871       768      402
  lonske:
  x     100      100    960300       560      282
  x     200      200   7721200       645      478
  x     300      300  16102700       630      393
*/
