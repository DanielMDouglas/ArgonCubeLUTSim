import sys, time, math
import argparse
import ROOT
import numpy as np


def LoadDataFile(path_data_file):

    print('\nLoadDataFile...')

    data = ROOT.TChain('events')
    data.Add(path_data_file)
    data.SetBranchStatus('*',0)
    data.SetBranchStatus('eventID',1)
    data.SetBranchStatus('event_nhits',1)
    data.SetBranchStatus('event_hits_x',1)
    data.SetBranchStatus('event_hits_y',1)
    data.SetBranchStatus('event_hits_z',1)
    data.SetBranchStatus('event_hits_q',1)

    print('using data file: ' + path_data_file.split('/')[-1] + ' (' + str(data.GetEntries()) + ' entries)')
    return data


def GetEntryListData(t_data, event):

    t_data.Draw('>>elist','eventID==%d' % event,'entrylist')

    return(ROOT.gDirectory.Get('elist'))


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


def GetVoxel(pos,lut_geometry):

    (lut_min,lut_max,lut_ndiv) = lut_geometry
    vox_xyz = np.floor(pos/(lut_max-lut_min)*lut_ndiv).astype(int)+lut_ndiv/2
    voxel = vox_xyz[2]*lut_ndiv[0]*lut_ndiv[1]+vox_xyz[1]*lut_ndiv[0]+vox_xyz[0]

    return voxel


def GetEntryListLUT(t_lut, voxel):

    t_lut.Draw('>>elist','Voxel==%d' % voxel,'entrylist')

    return(ROOT.gDirectory.Get('elist'))


def GetModuleCopy(pos):
    mod_x = (math.floor(pos[0]/(ModuleDimension[0]))+n_mod[0]/2)
    mod_z = (math.floor(pos[2]/(ModuleDimension[2]))+n_mod[2]/2)

    return int(mod_x+mod_z*n_mod[0])


def GetHalfDetCopy(pos):
    tpc_x = math.floor(pos[0]/(ModuleDimension[0]/2.))+n_mod[0]/2*n_tpc[0]

    return int(tpc_x)%2


def GetTPCActiveCenter(module_copynumber,halfDetector_copynumber):

    return TPCActiveCenter[str(module_copynumber)+str(halfDetector_copynumber)]


def GlobalToLUTCoord(pos,lut_geometry):

    # access LUT geometry
    (lut_min,lut_max,lut_ndiv) = lut_geometry

    # translate to TPCActive coordinates
    tpc_pos = pos-np.array(GetTPCActiveCenter(GetModuleCopy(pos),GetHalfDetCopy(pos)))

    # translate to TPCActive coordinates Module0
    tpc_pos = pos-np.array(GetTPCActiveCenter(GetModuleCopy(pos),GetHalfDetCopy(pos)))

    # negate drift and beam coordinates for mirrored TPCs
    if (GetHalfDetCopy(pos)==1):
        tpc_pos[0] = -tpc_pos[0]
        tpc_pos[2] = -tpc_pos[2]

    # translate to LUT coordinates
    lut_pos = tpc_pos+(lut_min+lut_max)/2

    return (lut_pos)


def TerminalPrintX(display,tphotons):
    op_channel_max = np.amax(tphotons)*2

    for op_channel in range(n_op_channel/2):
        line = (n_op_channel/2 - op_channel%(n_op_channel/2) -1)*2

        for i in range(int(round((tphotons[op_channel,0]+tphotons[op_channel+n_op_channel/2,0])/op_channel_max*10))):
            display[line][i] = '#'
        
        for i in range(int(round((tphotons[op_channel,1]+tphotons[op_channel+n_op_channel/2,1])/op_channel_max*10))):
            display[line][-i-1] = '#'


    print('\n')
    
    for i in range(48):
        print('\n'),
        for j in range(24+20):
            if (j==10 or j==22 or j==34):
                print('|'),
            if (display[i][j]=='*'):
                print('\x1b[1;36m' + display[i][j] + '\033[0;0m'),
            elif (display[i][j]=='#'):
                print('\x1b[1;33m' + display[i][j] + '\033[0;0m'),
            else:
                print(display[i][j]),

    print('\n---> x (drift) [%d photons/#]' % op_channel_max)


def TerminalPrintY(display,tphotons):
    op_channel_max = np.amax(tphotons)*2

    for op_channel in range(n_op_channel):
        line = (n_op_channel/2 - op_channel%(n_op_channel/2) -1)*2

        if (op_channel < 24):
            for i in range(int(round((tphotons[op_channel,0]+tphotons[op_channel,1])/op_channel_max*10))):
                display[line][i] = '#'
            
        else:
            for i in range(int(round((tphotons[op_channel,0]+tphotons[op_channel,1])/op_channel_max*10))):
                display[line][-i-1] = '#'


    print('\n')
    
    for i in range(48):
        print('\n'),
        for j in range(24+20):
            if (j==10 or j==34):
                print('|'),
            if (display[i][j]=='*'):
                print('\x1b[1;36m' + display[i][j] + '\033[0;0m'),
            elif (display[i][j]=='#'):
                print('\x1b[1;33m' + display[i][j] + '\033[0;0m'),
            else:
                print(display[i][j]),

    print('\n---> z (beam) [%d photons/#]' % op_channel_max)


def PhotonYield(Q):
    #Q = const_A_Birk/const_W_ion*(E/(1.+(const_k_Birk/(const_rho*const_drift))*E)) # [e]
    E = 1./(const_A_Birk/(Q*const_W_ion)-(const_k_Birk/(const_rho*const_drift))) # [MeV]
    Nq = E/const_W_ion
    return (Nq-Q)


def Simulate(t_data,t_lut,lut_geometry,event_number,n_evt,l_time,l_tphotons):
    
    # get event
    entry_list_event = GetEntryListData(t_data,event_number)
    if entry_list_event.GetEntry(0) == -1:
        print('\n\nskipping event no. %d because it looks broken...' % event_number)
        return
    print('\nsimulating event no. %d of %d (tree entry %d)...' % (event_number, n_evt, entry_list_event.GetEntry(0)))

    t_data.GetEntry(entry_list_event.GetEntry(0))
    n_hits = t_data.event_nhits
    print('%d hits found...' % (n_hits))
  
    # data container
    time = np.full((n_op_channel,np.prod(n_tpc)),20.)
    tphotons = np.zeros((n_op_channel,np.prod(n_tpc)))

    # event display
    display_x = np.empty((48,24+20),dtype='str')
    display_x[:] = ' '
    display_y = np.empty((48,24+20),dtype='str')
    display_y[:] = ' '

    # loop edep positions
    for hit in range(n_hits):

        # global position
        pos = np.array([t_data.event_hits_z[hit],t_data.event_hits_y[hit],t_data.event_hits_x[hit]])

        # event display
        #display_x[47-(int(pos[1]/TPCActiveDimension[1]*48+24))][(int(pos[0]/(ModuleDimension[2])*24+12+10))] = '*'
        #display_y[47-(int(pos[1]/TPCActiveDimension[1]*48+24))][(int(pos[2]/(ModuleDimension[2])*24+12+10))] = '*'
       
        # shift module0 coordinates 
        pos[0] = pos[0] - ModuleDimension[0]/2
        pos[1] = pos[1] - 218.236
        pos[2] = pos[2] - ModuleDimension[2]/2
        sys.stdout.write('\r    current position: ' + str(pos) + ('(%.1f %%)' % ((hit+1)/float(n_hits)*100)))
        sys.stdout.flush()

        if (GetModuleCopy(pos)!=0):
            continue

        # LUT position
        lut_pos = GlobalToLUTCoord(pos,lut_geometry)
        
        # voxel containing LUT position
        voxel = GetVoxel(lut_pos,lut_geometry)

        # voxel LUT entry-list
        entry_list_voxel = GetEntryListLUT(t_lut,voxel)

        # loop LUT entry-list
        for entry in range(entry_list_voxel.GetN()):
            t_lut.GetEntry(entry_list_voxel.GetEntry(entry))

            # tpc
            tpc = GetHalfDetCopy(pos)

            # optical channel
            op_channel = t_lut.OpChannel

            # swap upstream <-> downstream SiPMs for mirrored TPCs
            if (GetHalfDetCopy(pos)==1):
                op_channel = (op_channel+n_op_channel/2)%n_op_channel

            # make optical channels unique
            #op_channel = op_channel + TPCActiveID[str(GetModuleCopy(pos))+str(GetHalfDetCopy(pos))]*n_op_channel

            n_photons = t_lut.Visibility*PhotonYield(t_data.event_hits_q[hit]*245)#1E+3)

            if (t_lut.T1 < time[op_channel,tpc]):
                time[op_channel,tpc] = t_lut.T1
            
            if (op_channel % 12) > 5:
                n_photons *= norm_lcm_acl

            tphotons[op_channel,tpc] += n_photons

    # store visibility
    for tpc in range(2):
        for op_channel in range(n_op_channel):
            exec('l_time.ch%02d[%d] = %d' % (op_channel,tpc,int(time[op_channel][tpc])))
            exec('l_tphotons.ch%02d[%d] = %d' % (op_channel,tpc,int(tphotons[op_channel][tpc])))
            #if (op_channel == 0 and time[op_channel][tpc]<20.):
            #    print('')
            #    print(time[op_channel][tpc])
            #    print(tphotons[op_channel][tpc])

    #TerminalPrintX(display_x,tphotons)
    #TerminalPrintY(display_y,tphotons)


def main():

    event_number = args.event
    
    # link experimental data
    t_data = LoadDataFile(args.data)

    # get number of events
    t_data.GetEntry(t_data.GetEntries()-1)
    n_evt = t_data.eventID

    # link PhotonLibraryData
    t_lut = LoadLUT(args.lut)

    # access LUT geometry
    lut_geometry = GetLutGeometry(args.lut)

    # declare optical channel structure
    channel_struct_string = 'struct ChannelList {'

    for op_channel in range(n_op_channel):
        channel_struct_string = channel_struct_string + ('int ch%02d[%d]; ' % (op_channel,np.prod(n_mod)*np.prod(n_tpc)))

    channel_struct_string = channel_struct_string + '};'
    
    ROOT.gInterpreter.Declare(channel_struct_string)

    # declare branch variables
    l_eventID = ROOT.std.vector('Int_t')([0])
    l_time = ROOT.ChannelList()
    l_tphotons = ROOT.ChannelList()

    # initialize branch variables
    for tpc in range(np.prod(n_mod)*np.prod(n_tpc)):
        for op_channel in range(n_op_channel):

            exec('l_time.ch%02d[%d]=20' % (op_channel,tpc))
            exec('l_tphotons.ch%02d[%d]=0' % (op_channel,tpc))

    # create output tree
    f = ROOT.TFile(args.data.split('/')[-1][0:-5] + '_optSim2x2.root','RECREATE')
    tree = ROOT.TTree('t_out','t_out')

    # declare branch structures
    channel_branch_string = ('ch00[%d]/I' % (np.prod(n_mod)*np.prod(n_tpc)))

    for op_channel in range(n_op_channel):
        if op_channel == 0:
            continue
        channel_branch_string = channel_branch_string + (':ch%02d[%d]' % (op_channel,np.prod(n_mod)*np.prod(n_tpc)))

    # add branches
    tree.Branch('l_eventID', l_eventID)
    tree.Branch('l_time',l_time,channel_branch_string)
    tree.Branch('l_tphotons',l_tphotons,channel_branch_string)
    
    # simulate events
    if (event_number == -999):
        for event_number in range(120): #n_evt):
            Simulate(t_data,t_lut,lut_geometry,event_number,n_evt,l_time,l_tphotons)
            l_eventID[0] = event_number
            tree.Fill()
    else:
        Simulate(t_data,t_lut,lut_geometry,event_number,n_evt,l_time,l_tphotons)
        l_eventID[0] = event_number
        tree.Fill()

    tree.Write()


#===============================================>
########## GLOBAL VARIABLE DECLARATION ##########
#===============================================>

n_mod = [2,1,2]
n_tpc = [2,1,1]

ModuleDimension = [670.,2022.414,670.]

TPCActiveDimension = [302.723,1241.1,620.3]

# ArgonCube 2x2
TPCActiveCenter = {
        '00':[-487.949,-218.236,-335.],
        '01':[-182.051,-218.236,-335.],
        '10':[182.051,-218.236,-335.],
        '11':[487.949,-218.236,-335.],
        '20':[-487.949,-218.236,335.],
        '21':[-182.051,-218.236,335.],
        '30':[182.051,-218.236,335.],
        '31':[487.949,-218.236,335.],
        }

TPCActiveID = {
        '00':0,
        '01':1,
        '10':2,
        '11':3,
        '20':4,
        '21':5,
        '30':6,
        '31':7
        }

TPCActiveCopyNo = {value:key for key, value in TPCActiveID.items()}

n_op_channel = 48

norm_lcm_acl = 14.9

# NEST Electronic Recoil (ER) Model
# NEST liquid argon mean yields note, E. Kozlova & J. Mueller on behalf of the NEST collaboration

const_ArDensity = 1.3954    # [g/cm^3]
const_WArgon = 23.6         # [eV]
const_EField = 500          # [V/cm]
const_Nq = 51.9             # [quanta/keV]

const_ER_a = 32.988-552.988/(15.5578+const_EField/(-4.7+0.025115*math.exp(const_ArDensity/0.265360653)))**0.208889
const_ER_b = 2.01952+20.9/(1.105+(const_EField/0.4)**4.55)**7.502
const_ER_c = 0.642039*(1000/const_WArgon+6.5*(5-0.5/(const_EField/1047.408)**0.01851))
const_ER_d = 10.3842
const_ER_p1 = 1
const_ER_p2 = 10.304
const_ER_p3 = 24.3509
const_ER_p4 = 0.10535
const_ER_p5 = 0.7
const_ER_LET = -2.11259
const_ER_DokeBirks = 1052.264+(1.415935E+10-1652.264)/(-5+(const_EField/0.328038)**1.74654)

# ICARUS Birks' model
# https://doi.org/10.3390/instruments5010002

const_A_Birk    = 0.800     # Birk model parameter
const_k_Birk    = 0.0486    # Birk model parameter ((kV/cm)(g/cm^2)/MeV)
const_W_ion     = 23.6E-6   # ionization work function (MeV)
const_drift     = 0.5       # drift HV (kV/cm)
const_rho       = 1.396     # LAr density (g/cm^2)

#===============================================<


if __name__ == '__main__':
    print('')
    print('=======================================================================')
    print('====>  ArgonCube 2x2 Light Simulation using Photon Look-up-table  <====')
    print('=======================================================================')
    print('Patrick Koller, LHEP, University of Bern')
    print('Contact: patrick.koller@lhep.unibe.ch')

    parser = argparse.ArgumentParser(description='Define photon look-up-table, datafile and event to simulate.')
    parser.add_argument('data', help='path to data file')
    parser.add_argument('lut', help='path to photon look-up table')
    parser.add_argument('-e', '--event', type=int, help='number of single event to process', default='-999')
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(1)

    main()

