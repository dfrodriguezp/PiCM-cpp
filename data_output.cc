#include "functions.h"
#include "/usr/include/hdf5/serial/hdf5.h"
#include <string>

const int steps = parameters::steps;
const int N = parameters::N;

void writeSpace(std::string filename, const int steps, const int N, double data[N][steps*2])
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

	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);
	status = H5Fclose(file_id);
}

void writeMesh(std::string filename, const int gp, const int dr)
{
	hid_t	file_id, dataset_id, dataspace_id;
	hsize_t	dims[2];
	herr_t	status;
	int		RANK = 2;
	int 	data[gp*gp][2];

	file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	dims[0] = gp;
	dims[1] = gp;

	dataspace_id = H5Screate_simple(RANK, dims, NULL);
	dataset_id = H5Dcreate2(file_id, "mesh", H5T_STD_I32BE, dataspace_id,
							H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	for (int i = 0; i < gp*gp; i++)
	{
		for (int j = 0; j < 2; j++)
			data[i][j] = 
	}

	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);
	status = H5Fclose(file_id);
}
