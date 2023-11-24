'''
Class to evaluate the barreloid thalamus model on NetPyNE

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
'''
##################################################################################################################################################################################################################

class Evaluate():
    def getPopGids(sim):
        # --- Find pop names
        pop_gids={}
        for pop_name_ in sim.net.pops.keys():
            if '@' in pop_name_:    pop=pop_name_.split('__')[0].split('@')[0]
            else:                   pop=pop_name_.split('__')[0]
            # Dictionary of cell gids by pop
            if pop not in pop_gids.keys(): pop_gids.update({pop:[]})
            pop_gids[pop]+=sim.net.pops[pop_name_].cellGids
        return pop_gids
    
    def getPopConvergence(sim):
        # --- Pop parameters
        pop_gids  = Evaluate.getPopGids(sim)
        pop_names = list(pop_gids.keys())

        # --- Convergence to a population
        pop_convergence = {target_pop:{} for target_pop in pop_names}
        for target_pop in pop_names:
            for cell_gid in pop_gids[target_pop]:
                try:
                    for conn in sim.net.cells[cell_gid].conns:
                        if conn['preGid']=='NetStim': continue # ignore NetStim conns
                        if (conn['synMech']=='esyn') or 'gap' in conn['synMech']:continue # skipping electrical synapses | those should be accounted separetely
                        pre_pop = conn['label'].split('|')[1].split('@')[0]
                        if pre_pop not in pop_convergence[target_pop].keys():   pop_convergence[target_pop].update({pre_pop:1})
                        else:                                                   pop_convergence[target_pop][pre_pop]+=1
                except:print('cell gid out of range: ', cell_gid)
        return pop_convergence

    def getPopDivergence(sim):
        # --- Pop parameters
        pop_gids=Evaluate.getPopGids(sim)
        pop_names = list(pop_gids.keys())

        # --- Convergence to a population
        pop_divergence  = {source_pop:{} for source_pop in pop_names}
        for source_pop in pop_names:
            for cell_gid in pop_gids[source_pop]:
                for conn in sim.net.cells[cell_gid].conns:
                    if conn['preGid']=='NetStim': continue # ignore NetStim conns
                    if (conn['synMech']=='esyn') or 'gap' in conn['synMech']:continue # skipping electrical synapses | those should be accounted separetely
                    pre_pop = conn['label'].split('|')[1].split('@')[0]
                    post_pop= conn['label'].split('|')[2].split('@')[0]
                    if post_pop not in pop_divergence[pre_pop].keys():      pop_divergence[pre_pop].update({post_pop:1})
                    else:                                                   pop_divergence[pre_pop][post_pop]+=1
        return pop_divergence
    
    def getCellConvergence(sim):
        pop_gids        = Evaluate.getPopGids(sim)
        pop_convergence = Evaluate.getPopConvergence(sim)
        cell_convergence={}
        for target_pop in pop_convergence.keys():
            cell_convergence.update({target_pop:{}})
            for pre_pop in pop_convergence[target_pop].keys():
                cell_convergence[target_pop].update({pre_pop:pop_convergence[target_pop][pre_pop]/len(pop_gids[target_pop])})
        return cell_convergence

    def getCellDivergence(sim):
        pop_gids        = Evaluate.getPopGids(sim)
        pop_divergence  = Evaluate.getPopDivergence(sim)
        cell_divergence={}
        for pre_pop in pop_divergence.keys():
            cell_divergence.update({pre_pop:{}})
            for post_pop in pop_divergence[pre_pop].keys():
                cell_divergence[pre_pop].update({post_pop:pop_divergence[pre_pop][post_pop]/len(pop_gids[pre_pop])})
        return cell_divergence
    
    def getConvergenceRatio(sim):
        cell_convergence  = Evaluate.getCellConvergence(sim)
        convergence_ratio = {target_pop:{} for target_pop in cell_convergence.keys()}
        for target_pop in cell_convergence.keys():
            for source_pop in cell_convergence[target_pop].keys(): convergence_ratio[target_pop].update({source_pop:cell_convergence[target_pop][source_pop]/sum(cell_convergence[target_pop].values())})
        return convergence_ratio
    
    def getDivergenceRatio(sim):
        cell_divergence  = Evaluate.getCellDivergence(sim)
        divergence_ratio = {source_pop:{} for source_pop in cell_divergence.keys()}
        for source_pop in cell_divergence.keys():
            for target_pop in cell_divergence[source_pop].keys(): divergence_ratio[source_pop].update({target_pop:cell_divergence[source_pop][target_pop]/sum(cell_divergence[source_pop].values())})
        return divergence_ratio



    def runAll(sim):
        # --- Pop parameters
        pop_gids=Evaluate.getPopGids(sim)
        pop_names = list(pop_gids.keys())

        # --- Pop  Convergence
        pop_convergence = Evaluate.getPopConvergence(sim)
        # --- Pop  Divergence
        pop_divergence  = Evaluate.getPopDivergence(sim)
        # --- Cell Convergence
        cell_convergence = Evaluate.getCellConvergence(sim)
        # --- Cell Divergence
        cell_divergence  = Evaluate.getCellDivergence(sim)

        return pop_gids,pop_names,pop_convergence,pop_divergence,cell_convergence,cell_divergence