from neuron import h

def loadCell(hocName, MorphoName):

    h.load_file('import3d.hoc')
    h.load_file('stdrun.hoc')
    MorphologyPath = '/home/fernando/Documents/SCx_model/O1_data_physiology/morphologies/ascii/'
    gid = 1
    h.load_file('/home/fernando/Documents/SCx_model/O1_data_physiology/emodels_hoc/' + hocName + '.hoc')

    cell = getattr(h,  hocName)(gid,MorphologyPath,MorphoName)

    return cell
