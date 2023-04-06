from lutSim import *

def array_dtype(nbins_time):
    return np.dtype([('vis','f4'), ('t0','f4'),('time_dist','f4',(int(nbins_time),))])

def MakeNumpyArray(t_lut, lut_geometry, nbins_time):
    (lut_min, lut_max, lut_ndiv) = lut_geometry
    print(lut_min)
    print(lut_max)
    print(lut_ndiv)

    # output array has shape (nVoxX, nVoxY, nVoxZ, nOpChan, 2)
    # last dimension is 0 = visibility, 1 = t0
    outputLUTarr = np.zeros((int(lut_ndiv[2]),
                             int(lut_ndiv[1]),
                             int(lut_ndiv[0]),
                             n_op_channel),
                             dtype=array_dtype(nbins_time))

    for i in range(int(lut_ndiv[0])):
        for j in range(int(lut_ndiv[1])):
            for k in range(int(lut_ndiv[2])):

                voxelxyz = [i, j, k]

                voxel = WrapVoxelInds(voxelxyz, lut_geometry)
                
                print (voxelxyz, voxel)
                entry_list_voxel = GetEntryListLUT(t_lut, voxel)

                print(entry_list_voxel)

                for entry in range(entry_list_voxel.GetN()):
                    t_lut.GetEntry(entry_list_voxel.GetEntry(entry))

                    opChan = t_lut.OpChannel

                    outputLUTarr[k, j, i, opChan, 0] = t_lut.Visibility
                    outputLUTarr[k, j, i, opChan, 1] = t_lut.T1

    return outputLUTarr

def GetVoxelInds(voxel, lut_geometry):
    (lut_min,lut_max,lut_ndiv) = lut_geometry
    voxel_xyz = (
        voxel % int(lut_ndiv[0]),
        (voxel//int(lut_ndiv[0])) % int(lut_ndiv[1]),
        (voxel//(int(lut_ndiv[0])*int(lut_ndiv[1]))) % int(lut_ndiv[2])
        )

    return voxel_xyz

def MakeNumpyArrayFast(t_lut, lut_geometry, nbins_time):
    (lut_min, lut_max, lut_ndiv) = lut_geometry
    print(lut_min)
    print(lut_max)
    print(lut_ndiv)

    # output array has shape (nVoxX, nVoxY, nVoxZ, nOpChan, 2)
    # last dimension is 0 = visibility, 1 = t0
    outputLUTarr = np.zeros((int(lut_ndiv[2]),
                             int(lut_ndiv[1]),
                             int(lut_ndiv[0]),
                             n_op_channel),
                             dtype=array_dtype(nbins_time))

    n = t_lut.GetEntries()
    print(f'Voxels in tree: {n}')
    for entry in range(n):
        t_lut.GetEntry(entry)
        i,j,k = GetVoxelInds(t_lut.Voxel, lut_geometry)
        if entry < 10 or (n > 10 and entry > n - 10):
            print(entry, (i,j,k), ':', t_lut.Voxel, t_lut.OpChannel, '|', t_lut.Visibility, t_lut.T1)
        elif entry == 10:
            print('...')
        elif entry > 10 and (entry % 10000 == 0):
            print(f'{entry}/{n}', end='\r')

        opChan = t_lut.OpChannel
        outputLUTarr[k, j, i, opChan]['vis'] = t_lut.Visibility
        outputLUTarr[k, j, i, opChan]['t0'] = t_lut.T1
        outputLUTarr[k, j, i, opChan]['time_dist'] = np.frombuffer(t_lut.Time.fArray, dtype='f4', count=nbins_time+2)[1:-1] # skip under- and over-flow bins

    return outputLUTarr

def main():
    
    # link PhotonLibraryData
    t_lut = LoadLUT(args.lut)

    # lookup number of time bins
    t_lut.GetEntry(0)
    nbins_time = int(t_lut.Time.GetNbinsX())

    # access LUT geometry
    lut_geometry = GetLutGeometry(args.lut)

    arr = MakeNumpyArrayFast(t_lut, lut_geometry, nbins_time)

    #np.save(args.output, arr)
    np.savez_compressed(args.output, arr=arr)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Define photon look-up-table, datafile and event to simulate.')
    parser.add_argument('lut',
                        help='path to photon look-up table')
    parser.add_argument('-o',
                        '--output',
                        help='path to output numpy look-up table')
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(1)

    main()
