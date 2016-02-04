from SimPEG import *
import simpegDCIP as DC
import matplotlib.pyplot as plt


def getSrcList(nElecs, aSpacing, in2D=False, plotIt=False):

    elocs = np.arange(0,aSpacing*nElecs,aSpacing)
    elocs -= (nElecs*aSpacing - aSpacing)/2
    space = 1
    WENNER = np.zeros((0,),dtype=int)
    for ii in range(nElecs):
        for jj in range(nElecs):
            test = np.r_[jj,jj+space,jj+space*2,jj+space*3]
            if np.any(test >= nElecs):
                break
            WENNER = np.r_[WENNER, test]
        space += 1
    WENNER = WENNER.reshape((-1,4))


    if plotIt:
        for i, s in enumerate('rbkg'):
            plt.plot(elocs[WENNER[:,i]],s+'.')
        plt.show()

    # Create sources and receivers
    i = 0
    if in2D:
        getLoc = lambda ii, abmn: np.r_[elocs[WENNER[ii,abmn]],0]
    else:
        getLoc = lambda ii, abmn: np.r_[elocs[WENNER[ii,abmn]],0, 0]
    srcList = []
    for i in range(WENNER.shape[0]):
        rx = DC.RxDipole(getLoc(i,1),getLoc(i,2))
        src = DC.SrcDipole([rx], getLoc(i,0),getLoc(i,3))
        srcList += [src]

    return srcList



def example(aSpacing=2.5, nElecs=10, plotIt=False):

    surveySize = nElecs*aSpacing - aSpacing
    cs = surveySize/nElecs/4

    mesh = Mesh.TensorMesh([
            [(cs,10, -1.3),(cs,surveySize/cs),(cs,10, 1.3)],
            [(cs,3, -1.3),(cs,3,1.3)],
    #         [(cs,5, -1.3),(cs,10)]
        ],'CN')
    if plotIt:
        mesh.plotGrid(showIt=True)

    srcList = getSrcList(nElecs, aSpacing, in2D=True)
    survey = DC.SurveyDC(srcList)
    problem = DC.ProblemDC_CC(mesh)
    problem.pair(survey)

    return mesh, survey, problem
