"""
Parallel runners.

"""



from screws.frozen import FrozenOnly
import inspect
import os
from root.config import *
from tools.deprecated.serial_runners.INSTANCES.matrix3d_input_runner import Matrix3dInputRunner





def RunnerDataReader(filename):
    """
    To do such job, we do not need multiple cores. One master core can do all things as we basically just read a file.

    1). We will not be able to retrieve all information, but just obtain enough data for, like, visualization.
    Therefore, we are not able to restart a task with the data reader.

    2). Whenever we code a new Runner, we need to include it in this function.

    :param filename: The name of the file to be read.
    :return: A instance of the correct runner class.
    """
    if rAnk == mAster_rank:
        with open(filename, 'r') as f:
            # so the header of a runner file must be like: <Runner>-<A particular runner classname>-<......>
            contents = f.readlines()
            c0 = contents[0]
            assert c0[1:7] == 'Runner' and '>-<' in c0 and c0.split('>-<')[1] != '', \
                f" <Runner> : it is not a Runner file! Its header is {c0} which does not match the standard format: " \
                f"<Runner>-<...A particular runner classname...>-<...something whatever else...>"
            runner_name = c0.split('>-<')[1]
    else:
        runner_name = None

    runner_name = cOmm.bcast(runner_name, root=mAster_rank)

    if runner_name == 'Matrix3dInputRunner':
        # Actually we read it with ``ParallelMatrix3dInputRunner``. Since the parallel version is extended from a serial
        # version, so something is old. And the way of reading is also special.
        if rAnk == mAster_rank:
            SR = Matrix3dInputRunner.readfile(filename)
        else:
            SR = None
        RDO = ParallelMatrix3dInputRunner(SR)
    else:
        raise NotImplementedError(f"No class <{runner_name}> to read to data file: {filename}.")

    RDO.___lock_iterate___ = True # we lock the runner such that it can not ``iterate`` as it does not have a solver.

    return RDO # RDO stands for: `Runner` but with `Data` `Only`.





#---------------- template of all runners ------------------------------------------------------------------------------
class ParallelRunner(FrozenOnly):

    """A template for all parallel runners."""
    def __init__(self):
        """"""
        self.___lock_iterate___ = False

    def iterate(self, *args, **kwargs):
        if self.___lock_iterate___:
            raise Exception('This parallel runner is locked; it can not run ``iterate`` function. This is probably'
                            ' because it is read from a file. So it lack a solver.')
        else:
            return self.___iterate___(*args, **kwargs)

    @property
    def visualize(self):
        return self.___visualize___

    @property
    def results(self):
        return self.___results___


    # A parallel runner must have the following methods and properties _________________________________________________

    def ___iterate___(self, *args, **kwargs):
        return NotImplementedError()

    @property
    def ___visualize___(self):
        return NotImplementedError()

    @property
    def ___results___(self):
        return NotImplementedError()






#--------- A PARTICULAR PARALLEL RUNNER --------------------------------- BELOW ----------------------------------------
class ParallelMatrix3dInputRunner(ParallelRunner):
    """"""

    def __init__(self, solver):
        """"""
        super().__init__()

        self._solver_ = solver

        if rAnk == mAster_rank:
            if solver.__class__.__name__ == 'Matrix3dInputRunner':
                # this is only for RunnerDataReader. Never use this manually
                self._SR_ = solver
            else:
                self._solver_source_code_ = inspect.getsource(solver)
                self._solver_dir_ = os.path.abspath(inspect.getfile(solver))
                self._SR_ = Matrix3dInputRunner(solver,
                                                solver_source_code=self._solver_source_code_,
                                                solver_dir=self._solver_dir_)
        else:
            self._SR_ = None
            self._slave_visualize_ = ___SlaveParallelMatrix3dInputRunnerVisualize___(self)

        self._freeze_self_()


    def ___iterate___(self, I1, I2, I3, criterion='standard', writeto=None, **kwargs):
        """

        :param I1: The 1st matrix to be passed to the solver.
        :param I2: The 2nd matrix to be passed to the solver.
        :param I3: The list or tuple of the 3rd variable to be passed to the solver
        :param criterion: The criterion of how to parse the inputs `I1`, `I2`, `I3`. It can be one of
            'standard' :
                The 'standard' criterion stands for that: `i0` and `i1` are main
                variables and they do not change along `i2`. So we need `i0` and `i1`
                to be iterable. And each `i0[.]` or `i1[.]` need to be iterable.

                For example: `i0` and `i1`:
                    i0 = [[1, 2, 3],
                          [4, 5],
                          [6, 7, 8, 8]]
                    i1 = [[0.5, 0.1, 0.3],
                          [0, 2],
                          [3, 4, -2, -3]]
                    i2 = [0, 0.15, 0.2]
                Note that shape(i0[k]) must be equal to shape(i1[k]).

        :param writeto: The key words variable to be passed to the solver.
        :param kwargs: The key words variable to be passed to the solver. For a ParallelMatrix3dInputRunner, the same
            kwargs will be used for all runs. It is not possible to customize kwargs for each run.
        :return:
        """

        if sIze == 1:
            self._SR_.iterate(I1, I2, I3, criterion=criterion, writeto=writeto, saveto=False, **kwargs)

        else:
            if rAnk == mAster_rank:
                self._SR_.iterate(I1, I2, I3, criterion=criterion, writeto=writeto, saveto=False, **kwargs)
            else:
                I, J, K = cOmm.recv(source=mAster_rank, tag=rAnk+0.1) # position mark 1 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                cOmm.barrier()

                for k in range(K):  # we let the axis2 go at the last.
                    for i in range(I):  # we let the axis0 go secondly.
                        for j in range(J):  # we let the axis1 go firstly.

                            Compute_Or_Not = cOmm.recv(source=mAster_rank, tag=rAnk + 0.2)  # position mark 2 <<<<<<<<<<
                            cOmm.barrier()

                            if Compute_Or_Not:
                                INPUTS = cOmm.recv(source=mAster_rank, tag=rAnk + 0.3)  # position mark 3 <<<<<<<<<<<<<<
                                cOmm.barrier()
                                _ = self._solver_(INPUTS[0], INPUTS[1], INPUTS[2], **INPUTS[3])

                # we do nothing after all computation in slave cores.

    @property
    def ___visualize___(self):
        if rAnk == mAster_rank:
            return self._SR_.visualize
        else:
            return self._slave_visualize_

    @property
    def ___results___(self):
        if rAnk == mAster_rank:
            R = self._SR_.rdf
        else:
            R = None

        R = cOmm.bcast(R, root=mAster_rank)

        return R





# noinspection PyUnusedLocal
class ___SlaveParallelMatrix3dInputRunnerVisualize___(FrozenOnly):
    """We have this just to make that we can call visualize without indicate rAnk."""
    def __init__(self, pm3ir):
        self._pm3ir_ = pm3ir
        self._quick_ = ___SPM3IRV_quick___()
        self._freeze_self_()

    @property
    def quick(self):
        """ Access to the quick visualization methods."""
        return self._quick_

    @staticmethod
    def plot(*args, **kwargs):
        return None

    @staticmethod
    def semilogx(*args, **kwargs):
        return None

    @staticmethod
    def semilogy(*args, **kwargs):
        return None

    @staticmethod
    def loglog(*args, **kwargs):
        return None


# noinspection PyUnusedLocal
class ___SPM3IRV_quick___(FrozenOnly):
    def __init__(self):
        self._freeze_self_()

    @staticmethod
    def scatter(*args, **kwargs):
        return None


#--------- A PARTICULAR PARALLEL RUNNER --------------------------------- BELOW ----------------------------------------

#========= ALL PARTICULAR PARALLEL RUNNERS ==================== ABOVE ==================================================