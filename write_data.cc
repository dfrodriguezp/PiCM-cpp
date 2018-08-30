#include "parameters.h"

void writeData(std::string filename, std::string sname, const VecVec& data) {
    hid_t   file_id, dataset_id, dataspace_id, memspace_id, dcpl_id;
    hsize_t dims[2], dimsm[1], chunk_dims[2];
    hsize_t offset[2];
    hsize_t count[2];
    hsize_t stride[2];
    hsize_t block[2];
    Index   NX = data.size();
    Index   NY = data.at(0).size();
    Index   szip_options_mask, pixels_per_block;

    file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    dims[0] = NX;
    dims[1] = NY;

    dataspace_id = H5Screate_simple(2, dims, NULL);
    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    szip_options_mask = H5_SZIP_NN_OPTION_MASK;
    pixels_per_block = 10;
    H5Pset_szip(dcpl_id, szip_options_mask, pixels_per_block);  

    if (sname == "mesh" || sname == "energy") {
        chunk_dims[0] = Index(NX / 1);
        chunk_dims[1] = 1;
    } else {
        chunk_dims[0] = 1;
        chunk_dims[1] = Index(NY / 1);
    }

    H5Pset_chunk(dcpl_id, 2, chunk_dims);
    dataset_id = H5Dcreate2(file_id, sname.c_str(), H5T_IEEE_F64LE, dataspace_id,
                            H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    dimsm[0] = NY;
    memspace_id = H5Screate_simple(1, dimsm, NULL);

    count[0] = 1;
    count[1] = NY;

    stride[0] = 1;
    stride[1] = 1;

    block[0] = 1;
    block[1] = 1;

    for (Index i = 0; i < NX; i++) {
        offset[0] = i;
        offset[1] = 0;
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, stride, count, block);
        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, data.at(i).data());
    }

    // dataset_id = H5Dcreate2(file_id, sname.c_str(), H5T_IEEE_F64LE, dataspace_id,
    //                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // for (Index i = 0; i < NX; i++) {
    //     H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.at(i).data());
    // }    

    std::cout << sname << " written" << std::endl;
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);
    H5Pclose(dcpl_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
}