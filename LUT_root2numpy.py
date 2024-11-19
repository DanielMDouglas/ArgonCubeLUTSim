#import sys, time, math
import argparse
import ROOT
import numpy as np

def array_dtype(nbins_time):
    return np.dtype([('vis','f4'), ('t0','f4'), ('t0_avg','f4'), ('time_dist','f4',(int(nbins_time),))])

def LoadLUT(path_lut):

    print('\nLoadLUT...')

    lut = ROOT.TChain('PhotonLibraryData')
    lut.Add(path_lut)
    lut.SetBranchStatus('*',1)

    print('using lut: ' + path_lut.split('/')[-1] + ' (' + str(lut.GetEntries()) + ' entries)')

    return lut

def GetLutGeometry(path_lut):

    f = ROOT.TFile.Open(path_lut)
    lut_min = np.array([f.Min[0],f.Min[1],f.Min[2]])
    lut_max = np.array([f.Max[0],f.Max[1],f.Max[2]])
    lut_ndiv = np.array([f.NDivisions[0],f.NDivisions[1],f.NDivisions[2]])
    f.Close()

    return np.array([lut_min,lut_max,lut_ndiv])

def GetVoxelInds(voxel, lut_geometry):
    (lut_min,lut_max,lut_ndiv) = lut_geometry
    voxel_xyz = (
        voxel % int(lut_ndiv[0]),
        (voxel//int(lut_ndiv[0])) % int(lut_ndiv[1]),
        (voxel//(int(lut_ndiv[0])*int(lut_ndiv[1]))) % int(lut_ndiv[2])
        )

    return voxel_xyz

def MakeNumpyArrayFast(t_lut, lut_geometry, nbins_time, n_op_channel):
    (lut_min, lut_max, lut_ndiv) = lut_geometry
    print("LUT min: ", lut_min)
    print("LUT max: ", lut_max)
    print("LUT number of bins: ", lut_ndiv)

    # output array has shape (nVoxX, nVoxY, nVoxZ, nOpChan)
    outputLUTarr = np.zeros((int(lut_ndiv[2]),
                             int(lut_ndiv[1]),
                             int(lut_ndiv[0]),
                             int(n_op_channel)),
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
        time_dist = np.frombuffer(t_lut.Time.fArray, dtype='f4', count=nbins_time+2)[1:-1] # skip under- and over-flow bins
        if np.sum(time_dist) != 0:
            outputLUTarr[k, j, i, opChan]['time_dist'] = time_dist / np.sum(time_dist)
        else:
            outputLUTarr[k, j, i, opChan]['time_dist'] = time_dist
        time_dist_idx = np.arange(time_dist.shape[0])
        outputLUTarr[k, j, i, opChan]['t0_avg'] = np.sum(np.multiply(outputLUTarr[k, j, i, opChan]['time_dist'], time_dist_idx))

    return outputLUTarr

def main():
    
    # link PhotonLibraryData
    t_lut = LoadLUT(args.lut)

    # lookup number of time bins
    t_lut.GetEntry(0)
    nbins_time = int(t_lut.Time.GetNbinsX())
    print("number of time bins in the photon arrival distribution: ", nbins_time)

    # access LUT geometry
    lut_geometry = GetLutGeometry(args.lut)

    n_op_channel = t_lut.GetEntries()/np.prod(lut_geometry[2])
    print("number of optical channels: ", n_op_channel)

    arr = MakeNumpyArrayFast(t_lut, lut_geometry, nbins_time, n_op_channel)

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
