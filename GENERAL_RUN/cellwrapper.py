from neuron import h

def loadCell(hocName, MorphoName):

    h.load_file('import3d.hoc')
    h.load_file('stdrun.hoc')
    MorphologyPath = '/home/fernando/NetPyNE_SONATA/GENERAL_RUN/morphology/'
    gid = 1
    h.load_file('/home/fernando/NetPyNE_SONATA/SCx_model/O1_data_physiology/emodels_hoc/' + hocName + '.hoc')

    cell = getattr(h,  hocName)(gid,MorphologyPath,MorphoName)

    return cell
    
def loadCell2(hocName, MorphoName):

    h.load_file('import3d.hoc')
    h.load_file('stdrun.hoc')
    MorphologyPath = '/home/fernando/NetPyNE_SONATA/GENERAL_RUN/morphology/'
    gid = 1
    h.load_file('/home/fernando/NetPyNE_SONATA/GENERAL_RUN/cell.hoc')

    cell = getattr(h,  hocName)(gid,MorphologyPath+MorphoName)

    return cell
    
    
def loadCell3(hocName, MorphoName):    

    h.load_file("createsimulation.hoc")
    h.create_cell()
    cell = h.cell

    return cell
