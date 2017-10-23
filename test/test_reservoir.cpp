#include "random.h"

#include <algorithm>
#include <ctime>
#include <iostream>
#include <memory>


typedef unsigned int s_t;



double time_diff(clock_t t0, clock_t t1)
{
    return double(t1 - t0) / CLOCKS_PER_SEC;
}



void print_usage(std::string const & cmd, const double alpha, const unsigned s, const int v)
{
    std::cout
        << "usage: " << cmd << std::endl
        << "         --alpha alpha (default " << alpha << ")" << std::endl
        << "         --cap  capacity  (required)" << std::endl
        << "         -s  seed  (default " << s << ", for random)" << std::endl
        << "         -v  verbosity  (default " << v << ")" << std::endl;
}



int main(int argc, char ** argv)
{
    double alpha = 1.0;
    unsigned seed = 0;
    int verbose = 1;
    int capacity = 0;


    int iarg = 1;
    while (iarg < argc)
    {
        std::string arg{argv[iarg]};
        ++iarg;
        if (arg.compare("--alpha") == 0)
        {
            alpha = atof(argv[iarg]);
            assert(alpha >= 0.);
        } else if (arg.compare("--cap") == 0)
        {
            capacity = atoi(argv[iarg]);
            assert(capacity > 0);
        } else if (arg.compare("-s") == 0)
        {
            seed = atoi(argv[iarg]);
        } else if (arg.compare("-v") == 0)
        {
            verbose = atoi(argv[iarg]);
        } else
        {
            print_usage(argv[0], alpha, seed, verbose);
            return -1;
        }
        iarg++;
    }

    if (capacity < 1)
    {
        print_usage(argv[0], alpha, seed, verbose);
        return -1;
    }

    if (seed == 0)
    {
        seed = global_randomize();
    } else
    {
        global_seed(seed);
    }

    std::cout << "Random seed set to " << seed << std::endl;

    weighted_reservoir reservoir(capacity, alpha);


    std::cout << "Reservoir initiated with capacity " << capacity << std::endl;

    const int n_max = capacity * 5;

    clock_t t0, t1;
    double run_time;

    if (verbose > 0)
    {
        std::cout << std::endl;
    }


    for (int repeat = 0; repeat < 5; ++repeat)
    {
        s_t n_provided = pick_a_number(0.1, 1.0) * n_max;
        auto old_size = reservoir.size();

        t0 = clock();
        reservoir.keep_n_append(n_provided);
        t1 = clock();
        run_time = time_diff(t0, t1);

        if (verbose > 0)
        {
            std::cout << "Took " << run_time << " seconds to add "
                << n_provided << " data points to reservoir of current size "
                << old_size
                << std::endl;

            if (verbose > 1)
            {
                std::cout << "  n_provided: " << n_provided << ";  n_kept: "
                    << reservoir.n_kept() << ";  n_added: " << reservoir.n_appended() << ";  grand_total: "
                    << reservoir.grand_total() << std::endl;
                if (verbose > 2)
                {
                    auto idx_current = reservoir.idx_current();
                    auto n_current = reservoir.size();
                    std::cout << "  data points in reservoir: ";
                    for (s_t i = 0; i < n_current; ++i)
                    {
                        std::cout << " " << *idx_current++;
                    }
                    std::cout << std::endl;
                    if (verbose > 3)
                    {
                        std::cout << "      Points kept: ";
                        auto idx_kept = reservoir.idx_kept();
                        auto n_kept = reservoir.n_kept();
                        for (s_t i = 0; i < n_kept; ++i)
                        {
                            std::cout << " " << *idx_kept++;
                        }
                        std::cout << std::endl;
                        std::cout << "      Points appended: ";
                        auto idx_appended = reservoir.idx_appended();
                        auto n_appended = reservoir.n_appended();
                        for (s_t i = 0; i < n_appended; ++i)
                        {
                            std::cout << " " << *idx_appended++;
                        }
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;
                }
            }
        }
    }

    reservoir.clear();

    global_seed(seed);


    if (verbose > 0)
    {
        std::cout << std::endl;
    }



    for (int repeat = 0; repeat < 5; ++repeat)
    {
        s_t n_provided = pick_a_number(0.1, 1.) * n_max;

        auto old_size = reservoir.size();

        t0 = clock();
        reservoir.remove_n_inject(n_provided);
        t1 = clock();
        run_time = time_diff(t0, t1);

        if (verbose > 0)
        {
            std::cout << "Took " << run_time << " seconds to add "
                << n_provided << " data points to reservoir of current size "
                << old_size
                << std::endl;

            if (verbose > 1)
            {
                std::cout << "  n_provided: " << n_provided << ";  n_removed: "
                    << reservoir.n_removed() << ";  n_added: " << reservoir.n_injected() << ";  grand_total: "
                    << reservoir.grand_total() << std::endl;
                if (verbose > 2)
                {
                    auto idx_current = reservoir.idx_current();
                    auto n_current = reservoir.size();
                    std::cout << "  data points in reservoir: ";
                    for (s_t i = 0; i < n_current; ++i)
                    {
                        std::cout << " " << *idx_current++;
                    }
                    std::cout << std::endl;
                    if (verbose > 3)
                    {
                        std::cout << "      Points removed: ";
                        auto idx_removed = reservoir.idx_removed();
                        auto n_removed = reservoir.n_removed();
                        for (s_t i = 0; i < n_removed; ++i)
                        {
                            std::cout << " " << *idx_removed++;
                        }
                        std::cout << std::endl;
                        std::cout << "      Points injected: ";
                        auto idx_injected = reservoir.idx_injected();
                        auto n_injected = reservoir.n_injected();
                        for (s_t i = 0; i < n_injected; ++i)
                        {
                            std::cout << " " << *idx_injected++;
                        }
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;
                }
            }
        }
    }


    if (verbose > 0)
    {
        std::cout << std::endl;
    }


    t0 = clock();
    reservoir.export_to_file("reservoir.h5");
    t1 = clock();
    run_time = time_diff(t0, t1);

    if (verbose > 0)
    {
        std::cout << "Took " << run_time << " seconds to export reservoir of size "
            << reservoir.size() << " to file"
            << std::endl;
        if (verbose > 2)
        {
            auto n_current = reservoir.size();
            auto idx_current = reservoir.idx_current();
            std::cout << "  data points in reservoir: ";
            for (s_t i = 0; i < n_current; ++i)
            {
                std::cout << " " << *idx_current++;
            }
            std::cout << std::endl;
        }
    }



    weighted_reservoir reservoir_again;

    t0 = clock();
    reservoir_again.import_from_file("reservoir.h5");
    t1 = clock();
    run_time = time_diff(t0, t1);

    if (verbose > 0)
    {
        std::cout << "Took " << run_time << " seconds to import reservoir of size "
            << reservoir_again.size() << " from file"
            << std::endl;
        if (verbose > 2)
        {
            auto n_current = reservoir.size();
            auto idx_current = reservoir.idx_current();
            std::cout << "  data points in reservoir: ";
            for (s_t i = 0; i < n_current; ++i)
            {
                std::cout << " " << *idx_current++;
            }
            std::cout << std::endl;
        }
    }


    reservoir_again.export_to_file("reservoir_rewrite.h5");
        // Then use 'diff' to confirm the files 'reservoir.h5' and
        // 'reservoir_rewrite.h5' are identical.
}
