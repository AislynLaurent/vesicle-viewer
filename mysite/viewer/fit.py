class Fit:
    def __init__(self, project, sample):
        self.project = project
        self.sample = sample

        # code that is needed in any fit is set here
        self.zero_parameter = False
        self.data_exists = False
        self.datas 
    
    def get_fit_main(self):
        raise NotImplementedError("get_fit_main() not implemented")
    

class SymmetricalFit(Fit):
    def get_fit_main(self):
        return 0
    


class AsymmetricalFit(Fit):
    def get_fit_main(self):
        return 0

# generator function
def generate_fit_main(project, sample, *args, **kwargs):
    if project.model_type == "SM":
        fit = SymmetricalFit(project, sample, args, kwargs)
    elif project.model_type == "AS":
        fit = AsymmetricalFit(project, sample, args, kwargs)
    
    return fit.get_fit_main()