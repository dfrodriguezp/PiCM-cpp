#include "functions.h"
#include <string>

void writeSpace(std::string filename, const int N, const int steps, )
{
	hid_t	file_id, dataset_id, dataspace_id;
	hsize_t	dims[2];
	herr_t	status;
	int		RANK = 2;

	file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	dims[0] = N;
	dims[1] = steps * 2;

	dataspace_id = H5Screate_simple(RANK, dims, NULL);
	dataset_id = H5Dcreate2(file_id, "space", H5T_IEEE_F64BE, dataspace_id,
							H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}