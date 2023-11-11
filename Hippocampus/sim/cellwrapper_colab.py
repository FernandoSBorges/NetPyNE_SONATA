from neuron import h

def loadCell(hocName, MorphoName):

    h.load_file('import3d.hoc')
    h.load_file('stdrun.hoc')
    MorphologyPath = '/home/jovyan/work/NetPyNE-UI/CA1_netpyne/info/data-bbp/20191017/morphologies/swc/'
    gid = 1
    h.load_file('cells/hoc/' + hocName + '.hoc')

    cell = getattr(h,  hocName)(gid,MorphologyPath,MorphoName)

    return cell
