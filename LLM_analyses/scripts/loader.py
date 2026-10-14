import antimony
import roadrunner

def loada(path):
    antimony.clearPreviousLoads()
    code = antimony.loadAntimonyFile(path)
    if code < 0:
        raise RuntimeError(antimony.getLastError())
    module_name = antimony.getMainModuleName()
    sbml = antimony.getSBMLString(module_name)
    r = roadrunner.RoadRunner(sbml)
    return r
