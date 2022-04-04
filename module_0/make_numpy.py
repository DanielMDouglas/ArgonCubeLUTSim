from lutSim import *

def MakeNumpyArray(t_lut, lut_geometry):
    (lut_min, lut_max, lut_ndiv) = lut_geometry
    print(lut_min)
    print(lut_max)
    print(lut_ndiv)

    # output array has shape (nVoxX, nVoxY, nVoxZ, nOpChan, 2)
    # last dimension is 0 = visibility, 1 = t0
    outputLUTarr = np.zeros((int(lut_ndiv[2]),
                             int(lut_ndiv[1]),
                             int(lut_ndiv[0]),
                             n_op_channel, 2))

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

def main():
    
    # link PhotonLibraryData
    t_lut = LoadLUT(args.lut)

    # access LUT geometry
    lut_geometry = GetLutGeometry(args.lut)

    arr = MakeNumpyArray(t_lut, lut_geometry)

    np.save(args.output, arr)

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
