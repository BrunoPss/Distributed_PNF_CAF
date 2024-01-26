// C++ standard library includes
#include <algorithm>
#include <cstdio>
#include <ctime>
#include <chrono>
#include <time.h>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// CAF includes
#include "caf/all.hpp"
#include "caf/io/all.hpp"

// Boost includes
CAF_PUSH_WARNINGS
#ifdef CAF_GCC
#  pragma GCC diagnostic ignored "-Wdeprecated-copy"
#endif
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/random.hpp>
CAF_POP_WARNINGS

// Own includes
#include "int512_serialization.hpp"
#include "is_probable_prime.hpp"
#include "types.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::unordered_map;
using std::vector;

using boost::multiprecision::gcd;
using boost::multiprecision::int512_t;

using namespace caf;

namespace {

struct config : actor_system_config {
  string host = "localhost";
  uint16_t port = 0;
  size_t num_workers = 0;
  string mode;
  config() {
    opt_group{custom_options_, "global"}
      .add(host, "host,H", "server host (ignored in server mode)")
      .add(port, "port,p", "port")
      .add(num_workers, "num-workers,w", "number of workers (in worker mode)")
      .add(mode, "mode,m", "one of 'server', 'worker' or 'client'");
  }
};

// -- SERVER -------------------------------------------------------------------

void run_server(actor_system& sys, const config& cfg) {
  if (auto port = sys.middleman().publish_local_groups(cfg.port))
    cout << "published local groups at port " << *port << '\n';
  else
    cerr << "error: " << caf::to_string(port.error()) << '\n';
  cout << "press any key to exit" << std::endl;
  getc(stdin);
}

// -- CLIENT -------------------------------------------------------------------

// Client state, keep track of factors, time, etc.
struct client_state {
  // Joined group
  group grp;

  // Time count
  // CPU Time
  clock_t total_cpuTime = 0;
  // Real Time
  std::chrono::time_point<std::chrono::high_resolution_clock> real_time_start;
  std::chrono::time_point<std::chrono::high_resolution_clock> real_time_end;
  std::chrono::microseconds real_time;

  // Rho Cycle Counter
  int512_t rho_cycle_counter = 0;

  // Factor List
  vector<int512_t> factors;
  // Number to Factorize
  int512_t N = 0;
  // N included in the challenge packet
  int512_t sender_N = 0;
  int512_t quest_N = 0;
  // Factorization Fase
  int factorize_director = 0;

  // Getters
  int512_t get_questN() { return quest_N; }
  std::chrono::time_point<std::chrono::high_resolution_clock> get_realTimeStart() { return real_time_start; }
  std::chrono::time_point<std::chrono::high_resolution_clock> get_realTimeEnd() { return real_time_end; }
  std::chrono::microseconds get_realTime() { return real_time; }
  clock_t get_totalCpuTime() { return total_cpuTime; }
  int512_t get_rhoCycleCounter() { return rho_cycle_counter; }

  // Setters
  void set_questN(int512_t n) { quest_N = n; }
  void set_realTimeEnd(std::chrono::time_point<std::chrono::high_resolution_clock> end) { real_time_end = end; }
  void set_realTime(std::chrono::microseconds time) { real_time = time; }
  void add_cpuTime(clock_t t) { total_cpuTime += t; }
  void add_rhoCycleCounter(int512_t cycleCounter) { rho_cycle_counter += cycleCounter; }

  int512_t ask_number(stateful_actor<client_state>* self, caf::group grp) {
      // State Initialization
      factors.clear();
      quest_N=0;
      factorize_director=0;

      string str;
      do {
          cout << "Insert number to factorize" << endl;
          cout << "-> ";
          std::cin >> str;

          // Check exit command
          if (!str.compare("exit")) {
              self->leave(grp);
              exit(0);
          }
          // Check if contains only digits
          if (!std::all_of(str.begin(), str.end(), ::isdigit)) {
              cout << "[ERROR] -> Insert valid number!" << endl;
          }
          else {
              // Check if number is prime
              if (is_probable_prime(int512_t(str))) {
                  cout << "[INFO] -> Number " << int512_t(str) << " is prime!" << endl;
              }
          }
      } while (!std::all_of(str.begin(), str.end(), ::isdigit) || is_probable_prime(int512_t(str)));

      N = int512_t(str);
      sender_N = N;

      // Elapsed Time counter start (Wall Clock Time)
      real_time_start = std::chrono::high_resolution_clock::now();
      // Reset CPU Time count
      total_cpuTime = 0;
      // Reset Rho Cycle Count
      rho_cycle_counter = 0;

      return N;
  }

  // Manages the factorization process
  // Coordinates message sending to workers (tasks)
  void work_manager(stateful_actor<client_state>* self, int512_t num) {
      if (num == 1) {
          cout << "Number 1 is not prime and not factorisable" << endl;
          return;
      }

      // Checks if number is prime
      if (is_probable_prime(num)) {
          // Inserts the number in the solution list
          factors.push_back(num);
          //cout << "Push " << num << endl;
          factorize_director = 2;
          return;
      }
      if (factorize_director == 2) {
          sender_N = num;
      }
      factorize_director = 1;
      
      cout << "[REQUEST] --> " << num << endl;
      quest_N = num;
      self->send(self->state.grp, num);
  }
};

behavior client(stateful_actor<client_state>* self, caf::group grp) {
  // Join group and save it to send messages later.
  self->join(grp);
  self->state.grp = grp;
  
  cout << "Client Behavior" << endl;

  // Ask the user for a number
  self->state.N = self->state.ask_number(self, grp);
  self->state.sender_N=self->state.N;

  // Check initial number is prime
  if (is_probable_prime(self->state.N)) {
      cout << "[INFO] -> Number " << self->state.N << " is prime!" << endl;
  }
  else {
      // Send Work Request
      self->state.work_manager(self, self->state.N);
  }

  return {
    // Main worker response processor
    [=](actor_id id, const int512_t& factor, const int512_t& quest_N, clock_t cpu_t, int512_t rhoCycleCounter) {

    // Checks if the response is related to the current task
    if (quest_N == self->state.get_questN()) {
        cout << "[RESPONSE] -> ID: " << id << "  Factor: " << factor << endl;
        // CPU time handling
        self->state.add_cpuTime(cpu_t);

        // Rho Cycle Counter Handling
        self->state.add_rhoCycleCounter(rhoCycleCounter);

        // Factorization Phase 1
        if (self->state.factorize_director == 1) {
            self->state.work_manager(self, factor);
        }
        // Factorization Phase 2
        if (self->state.factorize_director == 2) {
            self->state.work_manager(self, self->state.sender_N / factor);
        }

        // Calculation of the solution list sum
        int512_t sum = 1;
        for(vector<int512_t>::iterator it = self->state.factors.begin(); it != self->state.factors.end(); ++it) {
            sum *= *it;
        }

        // Checks if the factorization is complete
        if (sum == self->state.N) {
            // Factorization Complete
            // Elapsed Time counter stop (Wall Clock Time)
            self->state.set_realTimeEnd(std::chrono::high_resolution_clock::now());
            self->state.set_realTime(std::chrono::duration_cast<std::chrono::microseconds>(self->state.get_realTimeEnd() - self->state.get_realTimeStart()));

            // Resets quest number
            self->state.set_questN(0);

            // Factor list formator
            std::ostringstream aux;
            std::copy(self->state.factors.begin(), self->state.factors.end(),
                      std::ostream_iterator<int512_t>(aux, ","));

            // Result Presentation
            cout << "-----------------------------------------------------------" << endl;
            cout << "| Results" << endl;
            cout << "|" << endl;
            cout << "| Complete Factorization -> " << aux.str() << endl;
            cout << "| CPU Time Spent         -> " << self->state.get_totalCpuTime() << " ticks ( " << double(self->state.get_totalCpuTime()) / double(CLOCKS_PER_SEC) << " s )" << endl;
            cout << "| Sum Rho Cycle Passes   -> " << self->state.get_rhoCycleCounter() << " cycle passes" << endl;
            cout << "| Elapsed Time           -> " << self->state.get_realTime().count() << " microseconds ( " << double(self->state.get_realTime().count()) / double(1000000) << " s )" << endl;
            cout << "-----------------------------------------------------------" << endl;

            // Restart Client
            self->state.work_manager(self, self->state.ask_number(self, grp));
        }
    }
    },
    [=](int512_t N) {},
    [=](const group_down_msg&) {
      cerr << "FATAL: Server is down" <<  endl;
    }
  };
}

void run_client(actor_system& sys, const config& cfg) {
  if (auto eg = sys.middleman().remote_group("vslab", cfg.host, cfg.port)) {
    auto grp = *eg;
    sys.spawn(client, grp);
  }
  // Server does not exist
  else {
    cerr << "error: " << caf::to_string(eg.error()) << '\n';
  }
}

// -- WORKER -------------------------------------------------------------------

// State specific to each worker.
struct worker_state {
  // Joined group.
  group grp;

  // Client Requester Address
  strong_actor_ptr client_requester = nullptr;

  // CPU Time Count
  clock_t cpu_time;

  // Rho Cycle pass Count
  int512_t rho_cycle_counter = 0;

  // Getters
  clock_t get_cpuTime() { return cpu_time; }
  strong_actor_ptr get_clientRequester() { return client_requester; }
  int512_t get_rhoCycleCounter() { return rho_cycle_counter; }

  // Setters
  void set_cpuTime(clock_t time) { cpu_time = time; }
  void set_clientRequester(strong_actor_ptr act) { client_requester = act; }
  void set_rhoCycleCounter(int512_t num) { rho_cycle_counter = num; }

  //Random Number Generator
  int512_t random_generator(int512_t min, int512_t max) {
      std::random_device rd;
      std::mt19937_64 generator(rd());
      return min + (generator() % max);
  }

  // Pollard's Rho Algorithm
  int512_t pollards_rho(int512_t N, int512_t a, stateful_actor<worker_state>* self, const actor& user) {
    if (N % 2 == 0) {
        rho_cycle_counter++;
        return 2;
    }

    int512_t x = random_generator(1, N);
    int512_t y = x;
    int512_t p = 1;
    int512_t d;

    // Rho Cycle
    bool exit_cond = true;
    while (p == 1 && exit_cond) {
        // Checks if mailbox is not empty
        if (self->mailbox().size()) {
            // Analyses next mailbox element
            // Client MSG
            if (self->peek_at_next_mailbox_element()->content().match_element<int512_t>(0)) {
                // Client Requester
                // -> Exit Rho Cycle (stop processing)
                if (self->peek_at_next_mailbox_element()->sender == client_requester) {
                    exit_cond = false;
                }
                // Other Client
                // -> Continue Rho Cycle
                else {
                    exit_cond = true;
                }
            }
            // Worker MSG
            else {
                // Worker from the group
                // -> Exit Rho Cycle (stop processing)
                if (self->peek_at_next_mailbox_element()->content().get_as<int512_t>(2) == N) {
                    exit_cond = false;
                }
                // Worker outside the group
                // -> Continue Rho Cycle
                else {
                    exit_cond = true;
                }
            }
        }

        // Factor Processing
        // If solution isn't found in 400000 cycle passes -> change random increment (a)
        if (rho_cycle_counter % 400000 == 0 && rho_cycle_counter != 0) {
            a = random_generator(1, N);
            cout << "[UPDATE] -> Increment (a): " << a << endl;
        }
        x = (x*x + a) % N;
        y = (y*y + a) % N;
        y = (y*y + a) % N;
        d = (y - x) % N;
        p = gcd(d, N);

        rho_cycle_counter++;
    }
    return p;
  }
};

behavior worker(stateful_actor<worker_state>* self, caf::group grp, const actor& user) {
  // Join group and save it to send messages later.
  self->join(grp);
  self->state.grp = grp;

  cout << "Worker Behavior" << endl;

  return {
    // Main Task Messages Handler
    [=](int512_t N) {
      // Start Timer (CPU Time)
      self->state.set_cpuTime(clock());

      // Register current sender
      self->state.set_clientRequester(self->current_sender());

      cout << "[TASK] -> N: " << N << endl;

      int512_t factor;
      // Generate Initial increment values
      int512_t a = self->state.random_generator(1, N);
      factor = self->state.pollards_rho(N, a, self, user);

      // Mailbox is not empty
      // -> Check next mailbox element
      if (self->mailbox().size()) {
          // Analyses next mailbox element
          // Client MSG
          if (self->peek_at_next_mailbox_element()->content().match_element<int512_t>(0)) {
              // Other Client
              // -> Send Answer
              if (self->peek_at_next_mailbox_element()->sender != self->state.get_clientRequester()) {
                  // Stop Timer (CPU Time)
                  self->state.set_cpuTime(clock() - self->state.get_cpuTime());

                  cout << "[SOLUTION] -> ID: " << self->id() << "  Factor: " << factor << endl;
                  // Group Message Send (task solution)
                  self->send(self->state.grp, self->id(), factor, N, self->state.get_cpuTime(), self->state.get_rhoCycleCounter());
              }
          }
          // Worker MSG
          else {
              // Worker Outside Group
              // -> Send Answer
              if (self->peek_at_next_mailbox_element()->content().get_as<int512_t>(2) != N) {
                  // Stop Timer (CPU Time)
                  self->state.set_cpuTime(clock() - self->state.get_cpuTime());

                  cout << "[SOLUTION] -> ID: " << self->id() << "  Factor: " << factor << endl;
                  // Group Message Send (task solution)
                  self->send(self->state.grp, self->id(), factor, N, self->state.get_cpuTime(), self->state.get_rhoCycleCounter());
              }
          }
      }
      // Mailbox is empty
      // -> Send Answer
      else {
          // Stop Timer (CPU Time)
          self->state.set_cpuTime(clock() - self->state.get_cpuTime());

          cout << "[SOLUTION] -> ID: " << self->id() << "  Factor: " << factor << endl;
          // Group Message Send (task solution)
          self->send(self->state.grp, self->id(), factor, N, self->state.get_cpuTime(), self->state.get_rhoCycleCounter());
      }

      // Reset Rho Cycle Counter
      self->state.set_rhoCycleCounter(0);
    },
    [&](actor_id id, const int512_t& factor, const int512_t& quest_N, clock_t cpu_t, int512_t rhoCycleCounter) {},
    [=](const group_down_msg&) {
      cerr << "FATAL: Server is down" <<  endl;
    }
  };
}

void run_worker(actor_system& sys, const config& cfg) {
  if (auto eg = sys.middleman().remote_group("vslab", cfg.host, cfg.port)) {
    auto grp = *eg;

    //Spawn Workers
    scoped_actor self{sys};
    for (int i=0; i < cfg.num_workers; i++) {
      sys.spawn(worker, grp, self);
    }
  } else {
    cerr << "error: " << caf::to_string(eg.error()) << '\n';
  }
  sys.await_all_actors_done();
}

// -- MAIN ---------------------------------------------------------------------

// dispatches to run_* function depending on selected mode
void caf_main(actor_system& sys, const config& cfg) {
  // Check serialization implementation. You can delete this.
  auto check_roundtrip = [&](int512_t a) {
    byte_buffer buf;
    binary_serializer sink{sys, buf};
    assert(sink.apply(a));
    binary_deserializer source{sys, buf};
    int512_t a_copy;
    assert(source.apply(a_copy));
    assert(a == a_copy);
  };
  check_roundtrip(1234912948123);
  check_roundtrip(-124);

  int512_t n = 1;
  for (int512_t i = 2; i <= 50; ++i)
    n *= i;
  check_roundtrip(n);
  n *= -1;
  check_roundtrip(n);

  // Dispatch to function based on mode.
  using map_t = unordered_map<string, void (*)(actor_system&, const config&)>;
  map_t modes{
    {"server", run_server},
    {"worker", run_worker},
    {"client", run_client},
  };
  auto i = modes.find(cfg.mode);
  if (i != modes.end())
    (i->second)(sys, cfg);
  else
    cerr << "*** invalid mode specified" << endl;
}

} // namespace

CAF_MAIN(io::middleman, id_block::vslab)
