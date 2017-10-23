#include "hdf5util.h"

#include <cassert>



template<typename T>
hid_t h5get_mem_type()
{
    assert(false);
        // Only use the specialized types.
}



template<>
hid_t h5get_mem_type<char>()
{
    return H5T_NATIVE_CHAR;
}

template<>
hid_t h5get_mem_type<unsigned char>()
{
    return H5T_NATIVE_UCHAR;
}


template<>
hid_t h5get_mem_type<short>()
{
    return H5T_NATIVE_SHORT;
}

template<>
hid_t h5get_mem_type<unsigned short>()
{
    return H5T_NATIVE_USHORT;
}


template<>
hid_t h5get_mem_type<int>()
{
    return H5T_NATIVE_INT;
}

template<>
hid_t h5get_mem_type<unsigned int>()
{
    return H5T_NATIVE_UINT;
}


template<>
hid_t h5get_mem_type<long>()
{
    return H5T_NATIVE_LONG;
}

template<>
hid_t h5get_mem_type<unsigned long>()
{
    return H5T_NATIVE_ULONG;
}


template<>
hid_t h5get_mem_type<float>()
{
    return H5T_NATIVE_FLOAT;
}


template<>
hid_t h5get_mem_type<double>()
{
    return H5T_NATIVE_DOUBLE;
}


hsize_t h5get_array_npoints(hid_t loc_id, const char * name)
{
    hid_t d = H5Dopen(loc_id, name, H5P_DEFAULT);
    if (d < 0)
    {
        return d;
    }

    hid_t s = H5Dget_space(d);
    if (s < 0)
    {
        H5Dclose(d);
        return s;
    }

    hsize_t size = H5Sget_simple_extent_npoints(s);
        // size < 0 if unsuccessful.

    H5Sclose(s);
    H5Dclose(d);
    return size;
}


