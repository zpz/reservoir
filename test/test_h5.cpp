#include "hdf5util.h"

#include <algorithm>
#include <iostream>
#include <memory>

typedef unsigned long val_t;


void print_usage(std::string const & cmd)
{
    std::cout
        << "usage: " << cmd << std::endl
        << "         -n  number of data values (required)" << std::endl;
}


int main(int argc, char ** argv)
{
    int n{0};

    int iarg = 1;
    while (iarg < argc)
    {
        std::string arg{argv[iarg]};
        ++iarg;
        if (arg.compare("-n") == 0)
        {
            n = atoi(argv[iarg]);
        } else
        {
            print_usage(argv[0]);
            return -1;
        }
        iarg++;
    }

    if (n < 1)
    {
        print_usage(argv[0]);
        return -1;
    }


    std::unique_ptr<val_t[]> data{new val_t[n]};
    std::iota(data.get(), data.get() + n, 0);

    hid_t file_id;
    hsize_t dims[1];

    file_id = H5Fcreate("make_read.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dims[0] = n;
    h5make_dataset_number(file_id, "data", 1, dims, data.get());
    H5Fclose(file_id);


    file_id = H5Fopen("make_read.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    h5read_dataset_number(file_id, "data", data.get());
    H5Fclose(file_id);


    file_id = H5Fcreate("make_read_rewrite.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dims[0] = n;
    h5make_dataset_number(file_id, "data", 1, dims, data.get());
    H5Fclose(file_id);

    return 0;
}

