"""This script contains class that calculates the dynamics metrics for the model systems. 

Returns:
        list: returns metric values of switching density, fold_repression, and transition_time in a list
"""

# Cell Density at the end of growth phase. 

class CalculateDynamicsMetrics():
    """This class contains functions that calculate the dynamics metrics of quorum sensing systems.
    """
    def __init__(self, gfp, time, cell_density) -> None:
        self.gfp = gfp
        self.time = time
        self.cell_density = cell_density
        self._switch_density = self.switch_density(self.gfp, self.cell_density)
        self._fold_repression = self.fold_repression(self.gfp)
        self._transition_time = self.transition_time(self.gfp, self.time)
        self._metrics = [self._switch_density, self._fold_repression, self._transition_time]

    def get_dynamic_metrics(self):
        "returns the dynamic metrics list"
        return self._metrics

    @staticmethod
    def switch_density(gfp, N):
        "Cell density at which the gfp drops down to 90% of max from max gfp observed"
        max_gfp = max(gfp)
        max_gfp_index = list(gfp).index(max_gfp)
        max_90_id= min(range(max_gfp_index,len(gfp)), key=lambda i: abs(gfp[i]- (max_gfp*0.9)))
        act_density = N[max_90_id]
        return act_density

    @staticmethod
    def fold_repression(gfp):
        "Calculates the ration of max gfp to min gfp"
        max_gfp = max(gfp)
        max_gfp_index = list(gfp).index(max_gfp)
        min_gfp = min(gfp[max_gfp_index:])
        fold  = max_gfp/min_gfp
        return fold

    @staticmethod
    def transition_time(gfp,t):
        "Calculates the time taken to reach 1.15*min_gfp from max gfp. "
        max_gfp = max(gfp)
        max_gfp_index = list(gfp).index(max_gfp)
        min_gfp = min(gfp[max_gfp_index:])  
        min_10_id= min(range(max_gfp_index,len(gfp)), key=lambda i: abs(gfp[i]- (min_gfp*1.15)))
        trans_time = t[min_10_id] - t[max_gfp_index]
        return trans_time