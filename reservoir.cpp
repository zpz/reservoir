#include "reservoir.h"

#include "hdf5.h"
#include "hdf5_hl.h"
#include "hdf5util.h"

#include <algorithm>
#include <memory>
#include <random>
#include <tuple>



/// Simple functions

std::default_random_engine & global_urng()
{
    static std::default_random_engine u{};
    return u;
}




void global_seed(unsigned int s)
{
    global_urng().seed(s);
}





unsigned int global_randomize()
{
    static std::random_device rd{};
    auto s = rd();
    global_seed(s);
    return s;
}




unsigned int global_randomise()
{
    return global_randomize();
}




int pick_a_number(int from, int thru)
{
    static std::uniform_int_distribution<int> d{};
    using param_t = decltype(d)::param_type;
    return d( global_urng(), param_t{from, thru} );
}




double pick_a_number(double from, double upto)
{
    static std::uniform_real_distribution<double> d{};
    using param_t = decltype(d)::param_type;
    return d( global_urng(), param_t{from, upto} );
}




//////////// functions for weighted_reservoir   ////////////////

using qquad_t = std::tuple<size_t, max_size_t, double, double>;
        // < index, time, u, priority >


weighted_reservoir::weighted_reservoir()
{
}


weighted_reservoir::weighted_reservoir(
        const size_t cap,
        const double alph)
{
    assert(alph >= 0.);
    assert(cap > 0);
    _alpha = alph;
    _capacity = cap;
    _chosen_times = std::unique_ptr<max_size_t[]>{new max_size_t[cap]()};
    _chosen_u = std::unique_ptr<double[]>{new double[cap]()};
        // The trailing '()' default-initializes the allocated memory;
        // otherwise 'valgrind' can give very puzzling memory error
        // messages related to 'export_to_file'.
        // These initializations also make the values upon export to and
        // import from disk files definite, which is a good thing.
    _idx_kept_or_removed = std::unique_ptr<size_t[]>{new size_t[cap]};
    _idx_appended_or_injected = std::unique_ptr<size_t[]>{new size_t[cap]};
        // The two above are not written out to files, hence
        // 'valgrind' won't complain about them.
}




// 'clear' empties the reservoir but does not change its
// '_alpha' and '_capacity' settings; neither does it free the allocated
// spaces for the internal arrays---note that the sizes of
// these spaces should stay consistent with '_capacity'.
void weighted_reservoir::clear()
{
    _current_size = 0;
    _grand_total = 0;
    _ref_L = 0;
    _kept_or_removed = 0;
    _n_kept_or_removed = 0;
    _n_appended_or_injected = 0;
}




bool weighted_reservoir::empty() const
{
    return
        _current_size == 0 &&
        _grand_total == 0 &&
        _ref_L == 0;
}




double weighted_reservoir::alpha() const
{
    return _alpha;
}




size_t weighted_reservoir::capacity() const
{
    return _capacity;
}





size_t weighted_reservoir::size() const
{
    return _current_size;
}



max_size_t weighted_reservoir::grand_total() const
{
    return _grand_total;
}






// Add new data points to the reservoir with bookkeeping
// for the new data points; no sampling is involved b/c
// the new total does not exceed the reservoir's capacity.
void direct_inject(
        max_size_t * const chosen_times,
        double * const chosen_u,
        size_t current_size,
        max_size_t grand_total,
        const size_t n_provided
        )
{
    std::uniform_real_distribution<double> urd{0.0, 1.0};
    auto urng = global_urng();

    for (size_t i = 0; i < n_provided; ++i)
    {
        chosen_times[current_size] = grand_total;
            // The first one gets index '0'.
        chosen_u[current_size] = urd(urng);
        ++current_size;
        ++grand_total;
    }
}




// Added new data to the reservoir with sampling, b/c
// the current size plus new data exceeds the reservoir's capacity.
// Upon return, the workspace 'quad' contains most useful info to be
// used for further processing.
// The reservoir's state is barely changed within this function;
// changes will be made after returning from this function.
void sample_inject(
        max_size_t const * const chosen_times,
        double const * const chosen_u,
        const size_t current_size,
        const max_size_t grand_total,
        const size_t n_provided,
        const size_t capacity,
        const double alpha,
        max_size_t & _ref_L,
        qquad_t * const quad,
            // Pre-allocated workspace, size should be at least
            //   current_size + n_provided
            // Upon return, its content is used for subsequent
            // processing.
        const size_t quad_len
        )
{
    assert(quad_len > capacity);

    std::uniform_real_distribution<double> urd{0.0, 1.0};
    auto urng = global_urng();

    if (current_size  > 0)
    {
        _ref_L = *std::min_element(chosen_times, chosen_times + current_size);
    }

    double factor = 1.0 / (grand_total - _ref_L + n_provided);

    for (size_t i = 0; i < current_size; ++i)
    {
        quad[i] = std::make_tuple(
                i,     // Index in existing data.
                chosen_times[i],  // Grand index in entire history.
                chosen_u[i],
                std::pow((chosen_times[i] - _ref_L) * factor, alpha) / chosen_u[i]
                );
            // FIXME: if _ref_L has not changed recently,
            // some speed improvement is possible here, b/c the pow does
            // not change except for a scaling.
    }


    max_size_t idx_grand = grand_total;
    max_size_t ref_diff = idx_grand - _ref_L;
    size_t idx_0 = current_size;
    size_t idx_new = 0;

    while (idx_new < n_provided)
    {
        size_t idx = idx_0;
        while (idx < quad_len)
        {
            auto u = urd(urng);
            quad[idx] = std::make_tuple(
                    idx_new,
                    idx_grand,
                    u,
                    std::pow(ref_diff * factor, alpha) / u
                    );
            ++idx;
            ++idx_new;
            if (idx_new == n_provided)
                break;
            ++idx_grand;
            ++ref_diff;
        }

        // Place the 'capacity' number of elements with the largest 'pow/u'
        // value at the front; these are the elements to stay in the
        // reservoir.
        std::nth_element(
                quad,
                quad + capacity,
                quad + idx,
                [](qquad_t const & x, qquad_t const & y)
                {  return std::get<3>(x) > std::get<3>(y); });

        idx_0 = capacity;
    }
}



// if (current_size + n_provided <= capacity)
// {
//
//    n_added = n_provided
//    n_kept = current_size
//    n_deleted = 0
//  =>
//    n_kept + n_added <= capacity
//    n_deleted + n_added = n_provided <= capacity
//
// } else
// {
//   if (n_provided < capacity)
//   {
//
//      0 <= capacity - current_size <= n_added <= n_provided
//      0 < capacity - n_provided <= n_kept <= current_size
//      0 <= n_deleted <= current_size + n_provided - capacity
//    =>
//      n_kept + n_added = capacity
//      n_deleted + n_added
//          <= (current_size + n_provided - capacity) + n_provided
//          <  current_size + n_provided
//          <  current_size + capacity
//          <  2 * capacity
//   } else
//   {
//      0 <= capacity - current_size <= n_added <= capacity
//      0 <= n_kept <= current_size
//      0 <= n_deleted <= current_size
//    =>
//      n_kept + n_added = capacity
//      n_deleted + n_added
//          <= current_size + capacity
//          <= 2 * capacity
//   }
// }
//
// Overall, the safe bounds are
//   n_added <= capacity
//   n_kept <= current_size <= capacity
//   n_deleted <= current_size <= capacity
//   n_kept + n_deleted = current_size <= capacity


void weighted_reservoir::keep_n_append(
        const size_t n_provided
        )
{
    assert(n_provided > 0);
    assert(_grand_total + n_provided > _grand_total);
        // Guard against overfow of 'max_size_t'.

    if (_current_size + n_provided <= _capacity)
    {
        direct_inject(
                _chosen_times.get(), _chosen_u.get(),
                _current_size, _grand_total,
                n_provided);

        _n_kept_or_removed = _current_size;
            // Number kept.
        std::iota(_idx_kept_or_removed.get(), _idx_kept_or_removed.get() + _current_size, 0);

        _n_appended_or_injected = n_provided;
            // Number appended.
        std::iota(_idx_appended_or_injected.get(), _idx_appended_or_injected.get() + n_provided, 0);

        _kept_or_removed = 1;
            // keep_n_append

        _current_size += n_provided;
        _grand_total += n_provided;
        return;
    }


    size_t buffer_size = std::min(_current_size + n_provided, _capacity + _capacity + _capacity);
    std::unique_ptr<qquad_t[]> buffer{new qquad_t[buffer_size]};
    auto workspace = buffer.get();

    sample_inject(
            _chosen_times.get(),
            _chosen_u.get(),
            _current_size,
            _grand_total,
            n_provided,
            _capacity,
            _alpha,
            _ref_L,   // by reference
            workspace,
            buffer_size);


    size_t nn;

    nn = 0;
    for (size_t i = 0; i < _capacity; ++i)
    {
        auto k = std::get<1>(workspace[i]);
        if (k < _grand_total)
        {
            _chosen_times[nn] = k;
            _chosen_u[nn] = std::get<2>(workspace[i]);
            _idx_kept_or_removed[nn] = std::get<0>(workspace[i]);
            ++nn;
        }
    }
    _n_kept_or_removed = nn;
        // Number kept.

    nn = 0;
    for (size_t i = 0, j = _n_kept_or_removed; i < _capacity; ++i)
    {
        auto k = std::get<1>(workspace[i]);
        if (k >= _grand_total)
        {
            _chosen_times[j] = k;
            _chosen_u[j] = std::get<2>(workspace[i]);
            _idx_appended_or_injected[nn] = std::get<0>(workspace[i]);
            ++nn;
            ++j;
        }
    }
    _n_appended_or_injected = nn;
        // Number appended.

    _kept_or_removed = 1;
        // keep_n_append

    _current_size = _capacity;
    _grand_total += n_provided;
}





void weighted_reservoir::remove_n_inject(
        const size_t n_provided
        )
{
    assert(n_provided > 0);
    assert(_grand_total + n_provided > _grand_total);
        // Guard against overfow of 'max_size_t'.

    if (_current_size + n_provided <= _capacity)
    {
        direct_inject(
                _chosen_times.get(), _chosen_u.get(),
                _current_size, _grand_total,
                n_provided);

        _n_kept_or_removed = 0;
            // Number removed.

        _n_appended_or_injected = n_provided;
            // Number injected.
        std::iota(_idx_appended_or_injected.get(), _idx_appended_or_injected.get() + _n_appended_or_injected, 0);

        _kept_or_removed = 2;
            // remove_n_inject

        _current_size += n_provided;
        _grand_total += n_provided;

        return;
    }


    size_t buffer_size = std::min(_current_size + n_provided, _capacity + _capacity + _capacity);
    std::unique_ptr<qquad_t[]> buffer{new qquad_t[buffer_size]};
    auto workspace = buffer.get();

    sample_inject(
            _chosen_times.get(),
            _chosen_u.get(),
            _current_size,
            _grand_total,
            n_provided,
            _capacity,
            _alpha,
            _ref_L,  // by reference
            workspace,
            buffer_size);


    std::fill_n(_chosen_times.get(), _capacity, _grand_total);
        // '_chosen_times' of all pre-existing data are smaller than
        // '_grand_total', hence filling it with '_grand_total'
        // serves as a flag after the '_chosen_times' of retained
        // pre-existing data are restored, that is, after the following
        // block, elements of '_chosen_times' with value '_grand_total'
        // are holes evacuated.

    size_t nn;

    nn = _current_size;
        // Number removed.
    for (size_t i = 0; i < _capacity; ++i)
    {
        auto k = std::get<1>(workspace[i]);
        if (k < _grand_total)
        {
            _chosen_times[std::get<0>(workspace[i])] = k;
                // _chosen_u does not need restoration b/c it was not
                // erased.
            --nn;
        }
    }
    _n_kept_or_removed = nn;


    for (size_t i = 0, j = 0; i < _current_size; ++i)
    {
        if (_chosen_times[i] == _grand_total)
        {
            _idx_kept_or_removed[j++] = i;
        }
    }


    nn = 0;
        // Number injected.
    size_t i = 0;
    size_t j = 0;
    while (j < _n_kept_or_removed)
    {
        while (std::get<1>(workspace[i]) < _grand_total)
        {
            ++i;
        }
        _chosen_times[_idx_kept_or_removed[j]] = std::get<1>(workspace[i]);
        _chosen_u[_idx_kept_or_removed[j]] = std::get<2>(workspace[i]);
        _idx_appended_or_injected[nn] = std::get<0>(workspace[i]);
        ++nn;
        ++j;
        ++i;
    }

    j = _current_size;
    while (i < _capacity)
    {
        auto k = std::get<1>(workspace[i]);
        if (k >= _grand_total)
        {
            _chosen_times[j] = k;
            _chosen_u[j] = std::get<2>(workspace[i]);
            _idx_appended_or_injected[nn] = std::get<0>(workspace[i]);
            ++nn;
            ++j;
        }
        ++i;
    }

    _n_appended_or_injected = nn;


    _kept_or_removed = 2;
        // remove_n_inject

    _current_size = _capacity;
    _grand_total += n_provided;
}





size_t weighted_reservoir::n_kept() const
{
    if (_kept_or_removed == 1)
        return _n_kept_or_removed;
    else
        return 0;
}




size_t const * weighted_reservoir::idx_kept() const
{
    if (_kept_or_removed == 1)
        return _idx_kept_or_removed.get();
    else
        return nullptr;
}




size_t weighted_reservoir::n_appended() const
{
    if (_kept_or_removed == 1)
        return _n_appended_or_injected;
    else
        return 0;
}




size_t const * weighted_reservoir::idx_appended() const
{
    if (_kept_or_removed == 1)
        return _idx_appended_or_injected.get();
    else
        return nullptr;
}





size_t weighted_reservoir::n_removed() const
{
    if (_kept_or_removed == 2)
        return _n_kept_or_removed;
    else
        return 0;
}




size_t const * weighted_reservoir::idx_removed() const
{
    if (_kept_or_removed == 2)
        return _idx_kept_or_removed.get();
    else
        return nullptr;
}




size_t weighted_reservoir::n_injected() const
{
    if (_kept_or_removed == 2)
        return _n_appended_or_injected;
    else
        return 0;
}




size_t const * weighted_reservoir::idx_injected() const
{
    if (_kept_or_removed == 2)
        return _idx_appended_or_injected.get();
    else
        return nullptr;
}





max_size_t const * weighted_reservoir::idx_current() const
{
    if (_current_size > 0)
        return _chosen_times.get();
    else
        return nullptr;
}




herr_t weighted_reservoir::export_to_file(hid_t loc_id) const
{
    assert(_capacity > 0);

    hsize_t dims[1];
    herr_t status;

    dims[0] = 1;

    status = h5make_dataset_number(loc_id, "alpha", 1, dims, &_alpha);
    if (status < 0)
        return status;

    status = h5make_dataset_number(loc_id, "capacity", 1, dims, &_capacity);
    if (status < 0)
        return status;

    status = h5make_dataset_number(loc_id, "current_size", 1, dims, &_current_size);
    if (status < 0)
        return status;

    status = h5make_dataset_number(loc_id, "grand_total", 1, dims, &_grand_total);
    if (status < 0)
        return status;

    status = h5make_dataset_number(loc_id, "ref_L", 1, dims, &_ref_L);
    if (status < 0)
        return status;


    if (_capacity > 0)
    {
        dims[0] = _capacity;
        status = h5make_dataset_number(loc_id, "chosen_times", 1, dims, _chosen_times.get());
        if (status < 0)
            return status;
        status = h5make_dataset_number(loc_id, "chosen_u", 1, dims, _chosen_u.get());
        if (status < 0)
            return status;
    }


    return 0;
}




herr_t weighted_reservoir::export_to_file(hid_t loc_id, char const * name) const
{
    if (name[0] == '.' && name[1] == '\0')
    {
        return this->export_to_file(loc_id);
    } else
    {
        hid_t group_id = H5Gcreate(loc_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (group_id < 0)
        {
            return group_id;
        }
        auto status = this->export_to_file(group_id);
        H5Gclose(group_id);
        return status;
    }
}





herr_t weighted_reservoir::export_to_file(char const * file) const
{
    hid_t file_id = H5Fcreate(file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0)
    {
        return file_id;
    }
    herr_t status = this->export_to_file(file_id);
    H5Fclose(file_id);
    return status;
}




herr_t weighted_reservoir::import_from_file(hid_t loc_id)
{
    assert(this->empty());

    herr_t status;


    auto old_capacity = _capacity;

    status = h5read_dataset_number(loc_id, "alpha", &_alpha);
    if (status < 0)
        return status;

    status = h5read_dataset_number(loc_id, "capacity", &_capacity);
    if (status < 0)
        return status;

    status = h5read_dataset_number(loc_id, "current_size", &_current_size);
    if (status < 0)
        return status;

    status = h5read_dataset_number(loc_id, "grand_total", &_grand_total);
    if (status < 0)
        return status;

    status = h5read_dataset_number(loc_id, "ref_L", &_ref_L);
    if (status < 0)
        return status;



    assert(_capacity > 0);
        // This is guaranteed by 'export_to_file'.

    if (old_capacity != _capacity)
    {
        if (_chosen_times != nullptr)
        {
            _chosen_times.reset(nullptr);
        }
        _chosen_times = std::unique_ptr<max_size_t[]>{new max_size_t[_capacity]()};
        if (_chosen_u != nullptr)
        {
            _chosen_u.reset(nullptr);
        }
        _chosen_u = std::unique_ptr<double[]>{new double[_capacity]()};
    }

    status = h5read_dataset_number(loc_id, "chosen_times", _chosen_times.get());
    if (status < 0)
        return status;

    status = h5read_dataset_number(loc_id, "chosen_u", _chosen_u.get());
    if (status < 0)
        return status;



    _kept_or_removed = 0;
    _n_kept_or_removed = 0;
    _n_appended_or_injected = 0;
    if (old_capacity != _capacity)
    {
        if (_idx_kept_or_removed != nullptr)
        {
            _idx_kept_or_removed.reset(nullptr);
        }
        _idx_kept_or_removed = std::unique_ptr<size_t[]>{new size_t[_capacity]};
        if (_idx_appended_or_injected != nullptr)
        {
            _idx_appended_or_injected.reset(nullptr);
        }
        _idx_appended_or_injected = std::unique_ptr<size_t[]>{new size_t[_capacity]};
    }

    return 0;
}




herr_t weighted_reservoir::import_from_file(hid_t loc_id, char const * name)
{
    if (name[0] == '.' && name[1] == '\0')
    {
        return this->import_from_file(loc_id);
    } else
    {
        hid_t group_id = H5Gopen(loc_id, name, H5P_DEFAULT);
        if (group_id < 0)
        {
            return group_id;
        }
        auto status = this->import_from_file(group_id);
        H5Gclose(group_id);
        return status;
    }
}



herr_t weighted_reservoir::import_from_file(char const * file)
{
    assert(this->empty());
    hid_t file_id = H5Fopen(file, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0)
    {
        return file_id;
    }
    herr_t status = this->import_from_file(file_id);
    H5Fclose(file_id);
    return status;
}


