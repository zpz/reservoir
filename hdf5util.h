#ifndef HDF5UTIL_H
#define HDF5UTIL_H


// FIXME:
// Use version 2 of H5Dopen H5Acreate and H5Dcreate
// This makeshift solution is probably related to the fact
// that HDF5 version (1.8.4?) on Colossus is slightly lower
// than on my laptop, which is 1.8.9.
//
// NOTE: this should come BEFORE #include "hdf5.h".

#define H5Dopen_vers 2
#define H5Gopen_vers 2
#define H5Gcreate_vers 2
#define H5Dcreate_vers 2

#include "hdf5.h"
#include "hdf5_hl.h"

#include <cassert>
#include <type_traits>




template<typename T>
hid_t h5get_mem_type();



template<typename T>
hid_t h5get_disk_type()
{
    auto size = sizeof(T);
    if (std::is_floating_point<T>::value)
    {
        if (size == 4)
            return H5T_IEEE_F32LE;
        else if (size == 8)
            return H5T_IEEE_F64LE;
        else
            assert(false);
    } else if (std::is_integral<T>::value)
    {
        if (std::is_signed<T>::value)
        {
            if (size == 1)
                return H5T_STD_I8LE;
            else if (size == 2)
                return H5T_STD_I16LE;
            else if (size == 4)
                return H5T_STD_I32LE;
            else if (size == 8)
                return H5T_STD_I64LE;
            else
                assert(false);
        } else
        {
            if (size == 1)
                return H5T_STD_U8LE;
            else if (size == 2)
                return H5T_STD_U16LE;
            else if (size == 4)
                return H5T_STD_U32LE;
            else if (size == 8)
                return H5T_STD_U64LE;
            else
                assert(false);
        }
    } else
    {
        assert(false);
    }
}




template<typename T>
herr_t h5make_dataset_number(hid_t loc_id, const char * dset_name, int rank, const hsize_t *dims, const T * buffer)
{
    return H5LTmake_dataset(loc_id, dset_name, rank, dims,
            h5get_mem_type<T>(),
            static_cast<const void *>(buffer));
}



template<typename T>
herr_t h5read_dataset_number(hid_t loc_id, const char * dset_name, T * buffer)
{
    return H5LTread_dataset(loc_id, dset_name,
            h5get_disk_type<T>(),
            static_cast<void *>(buffer));
}


// Get the size (number of data entries) in a dataset that is known to
// be a 1D array.
hsize_t h5get_array_npoints(hid_t loc_id, const char * name);

#endif   // HDF5UTIL_H

