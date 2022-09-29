# -*- coding: utf-8 -*-
import os
import pandas as pd
from tqdm import tqdm
from time import time, sleep
from screws.freeze.main import FrozenClass
from screws.miscellaneous.timer import NumpyStyleDocstringReader
from screws.miscellaneous.timer import MyTimer, randomStringDigits
import inspect, pickle, psutil

from root.config.main import cOmm, rAnk, mAster_rank, ASSEMBLE_COST

from tools.iterators.base.monitor.main import IteratorMonitor
from tools.iterators.base.visualize.main import IteratorVisualize

class Iterator(FrozenClass):
    """A parent of all iterators.

    :param auto_save_frequency:
        After `auto_save_frequency` iterations, we will do the auto save to `RDF_filename`.csv
        And if auto_save_frequency = True, this frequency will be decided by the program according
        to `monitor_factor`.
    :param float monitor_factor:
        A float between 0 or >0.1. When it is 0, no monitor. When it is 1, normal frequency.
        When it > 1, higher frequency.
    :param RDF_filename:
        We will save the RDF to this file (as well as the auto save). And if it already exists, we
        first read RDF from it.
    :param name:
        The name of this iterator. And we also name the graph report picture file as:
        'MPI_IGR_' + name.
    :param save_to_mitr:
        Do we save the results into the `.mitr` file?
    """
    def __init__(self,
        auto_save_frequency = True,
        monitor_factor = 1,
        RDF_filename = None,
        real_time_monitor = False,
        name = None,
        save_to_mitr = False
        ):

        # these four will be initialized in the particular iterator.
        assert hasattr(self, '_t0_')
        assert hasattr(self, '_dt_')
        assert hasattr(self, '_max_steps_')
        assert hasattr(self, '_max_time_')
        # ... not needed to be initialized in the particular iterator.
        self._running_step_ = 1
        self._computed_steps_ = 0
        if rAnk == mAster_rank:
            self._monitor_ = IteratorMonitor(self,
                                             auto_save_frequency,
                                             RDF_filename,
                                             monitor_factor,
                                             real_time_monitor)
            if name is None:
                name ='Iterator-' + randomStringDigits(8) + '-' + str(id(self))[-5:]
        else:
            self._monitor_ = None
            name = None

        self._visualize_ = IteratorVisualize(self)

        name = cOmm.bcast(name, root=mAster_rank)
        self.standard_properties.name = name
        self._save_to_mitr_ = save_to_mitr
        self.___cpu_load___ = None
        self.___assembling_cost___ = None
        self._freeze_self_()

    @property
    def t0(self):
        # noinspection PyUnresolvedReferences
        return self._t0_

    @property
    def dt(self):
        """(float) Current ``dt``. So the next time will be ``self.t + self.dt``."""
        return self._dt_
    @dt.setter
    def dt(self,dt):
        assert dt > 0
        self._dt_ = dt

    @property
    def max_steps(self):
        # noinspection PyUnresolvedReferences
        return self._max_steps_
    @property
    def max_time(self):
        # noinspection PyUnresolvedReferences
        return self._max_time_

    @property
    def t(self):
        """(float) Current time."""
        return self._t_
    @t.setter
    def t(self, t):
        assert t > self.t
        self._t_ = t

    @property
    def running_step(self):
        return self._running_step_
    @running_step.setter
    def running_step(self, running_step):
        assert running_step == self.running_step + 1
        self._running_step_ = running_step

    @property
    def computed_steps(self):
        return self._computed_steps_
    @computed_steps.setter
    def computed_steps(self, computed_steps):
        assert computed_steps == self.computed_steps + 1
        self._computed_steps_ = computed_steps

    @property
    def monitor(self):
        """The monitor of this iterator."""
        return self._monitor_

    @property
    def visualize(self):
        return self._visualize_

    def __call__(self, solver, initial):
        """
        By calling the iterator, we apply the iterator to a particular solver (normally a function) and initialize
        the result DataFrame.

        :param solver:
        :param initial: The initial values for t0. will be put in the index 0 of the result DataFrame.
        :return:
        """
        ds = NumpyStyleDocstringReader(solver)
        solver_par = ds.Parameters
        solver_ret = ds.Returns

        assert len(solver_par) == 2 and solver_par[0] == 'tk' and solver_par[1] == 'tk1', \
            f" <iterator> : to use iterator, two parameters need to be 'tk', 'tk1'. {solver_par}"
        assert solver_ret[0] == 'exit_code', "First output must be exit_code."
        assert solver_ret[1] == 'shut_down', "Second output must be shut_down."
        assert solver_ret[2] == 'message', "Third output must be message."
        assert len(solver_ret) > 3, "need outputs."

        assert len(initial) == len(solver_ret) - 3, \
            f" len(initial) = {len(initial)} != (len(solver_ret) - 3)={len(solver_ret) - 3}"

        self._melt_self_()
        self._exit_code_ = 1
        self._shut_down_ = False
        self.___shut_down_reason___ = None
        self._message_ = None
        self._t_ = self.t0 # initialing `t` here

        self._solver_ = solver
        self._solver_source_code_ = inspect.getsource(solver)
        self._solver_dir_ = os.path.abspath(inspect.getfile(solver))

        if rAnk == mAster_rank:
            rdf_ret = solver_ret[3:]
            RDF = dict()
            RDF['t'] = self.t0
            RDF['dt'] = self.dt
            for k, key in enumerate(rdf_ret):
                RDF[key] = initial[k]

            if self.monitor.RDF_filename is None:
                self._RDF_ = pd.DataFrame(RDF, index=[0, ])
            else: # to restart a quest from a RDF file.
                try:
                    self.___PRIVATE_read_RDF___(self.monitor.RDF_filename)
                except NotImplementedError:
                    raise Exception(f' file: {self.monitor.RDF_filename} extension wrong!')
                except FileNotFoundError:
                    self._RDF_ = pd.DataFrame(RDF,index=[0,])
        else:
            self._RDF_ = None
        self._freeze_self_()

    def __next__(self):
        """Must be implemented in tha child classes."""
        raise NotImplementedError()

    @property
    def RDF(self):
        """(pandas.DataFrame) Result DataFrame."""
        return self._RDF_

    def ___PRIVATE_read_RDF___(self, fileName):
        """
        We renew the ``RDF`` from a file. We accept .pkl, .csv and .xlsx file.

        The iterations in this ``RDF`` will not be repeated. So this is kind of to resume a simulation.
        """
        FORMAT = fileName.split('.')[-1]
        if FORMAT == 'csv':
            self._RDF_ = pd.read_csv(fileName, index_col=0)
        else:
            raise NotImplementedError(f'Can only read .csv file, now get .{FORMAT} file')




    @property
    def shut_down(self):
        """(bool) Return ``True`` if the iterator is shut_down.

        ``shut down`` property is an interface for a particular solver to decide to stop or not when its max iterations
        are not reached yet.
        """
        assert self._shut_down_ in (1, 0, True, False), "Now, shut_down must be 1, 0, True or False."
        if self._shut_down_:
            assert isinstance(self.___shut_down_reason___, str), "Shut down reason need to be str, now is not!"
        else:
            assert self.___shut_down_reason___ is None
        return self._shut_down_

    @shut_down.setter
    def shut_down(self, shut_down):
        assert isinstance(shut_down, str) or shut_down in (1, 0, True, False), \
            "Shut down must be 1, 0, True, False or str. When str, put str to shut_down_reason, make shut_down = 1"

        if isinstance(shut_down, str):
            self.___shut_down_reason___ = shut_down
            self._shut_down_ = 1
        else:
            self._shut_down_ = shut_down
            if shut_down:
                self.___shut_down_reason___ = ''

        assert shut_down in (1, 0, True, False), "Now, shut_down must be 1, 0, True or False."

        if shut_down:
            assert isinstance(self.___shut_down_reason___, str), "Shut down reason need to be str, now is not!"
        else:
            assert self.___shut_down_reason___ is None

    @property
    def shut_down_reason(self):
        if isinstance(self.___shut_down_reason___, str):
            assert self.shut_down, "Must shut down to have a reason."
        else:
            assert self.___shut_down_reason___ is None
        return self.___shut_down_reason___




    @property
    def exit_code(self):
        """Return the exit code of the solver for the last run."""
        return self._exit_code_

    @exit_code.setter
    def exit_code(self, exit_code):
        self._exit_code_ = exit_code



    @property
    def message(self):
        """List(str) Return the messages of the solver for the last run."""
        return self._message_

    @message.setter
    def message(self, message):
        if isinstance(message, str):
            message = [message,]
        assert isinstance(message, (list, tuple)), "message must be str or list or tuple."
        for i, mi in enumerate(message):
            assert isinstance(mi, str), f"message must be tuple or list of str, " \
                                        f"now message[{i}]={mi} is not str"
        self._message_ = message


    def ___PRIVATE_append_outputs_to_RDF___(self, outputs):
        """"""
        self.RDF.loc[self.running_step-1] = [self.t, self.dt] + list(outputs[3:])



    #-------------- core methods: run and read ------------------------------------------------


    def run(self):
        """To run the iterator.

        To make it work properly, we have to make sure the solver return exactly the same
        outputs in all cores. Otherwise, cores may stop after different time steps. This very
        cause some serious problems.
        """
        if rAnk == mAster_rank:
            sleep(0.05)
            self.monitor._ft_firstRun_ = time()
            self.monitor._preparation_time_ = time() - self.monitor._ft_start_time_
            pbar = tqdm(total=self.max_steps,
                        desc = MyTimer.current_time() + ' <' + self.__class__.__name__ + '>')

        IN = 0 # if in the while loop.

        while not self.shut_down and self.running_step <= self.max_steps:

            # we now make sure all cores do together or no one does.

            IN = 1

            ALL_IN = cOmm.gather(IN, root=mAster_rank)

            # decide do the iteration or pass it (when running_step is in the RDF) ...
            if rAnk == mAster_rank:

                assert all(ALL_IN), "Not all cores are in the while loop."

                if self.running_step in self.RDF.index:
                    _pass = 1
                else:
                    _pass = 0
            else:
                _pass = None
            _pass = cOmm.bcast(_pass, root=mAster_rank)
            # Do or pass ...
            dt = self.dt # must do this since dt will be updated in next()

            if rAnk == mAster_rank:
                psutil.cpu_percent(None)
                self.___assembling_cost___ = 0

            if _pass:
                outputs = [1, 0, 'existing solution']
            else:
                outputs = next(self)

            if rAnk == mAster_rank:
                self.___cpu_load___ = psutil.cpu_percent(None)
                if len(ASSEMBLE_COST['recent']) > 0:
                    self.___assembling_cost___ = sum(ASSEMBLE_COST['recent'])
                    ASSEMBLE_COST['accumulated'] += self.___assembling_cost___
                    ASSEMBLE_COST['recent'] = list()

            # updates other self properties
            self.running_step += 1
            self.computed_steps += 1
            self.t += dt
            self.exit_code, self.shut_down, self.message = outputs[:3]
            if rAnk == mAster_rank:
                # update RDF ...
                if _pass == 0: # not already in RDF
                    self.___PRIVATE_append_outputs_to_RDF___(outputs)

                # update monitor ...
                self.monitor.do.update() # always update monitor ...

                if _pass == 0: # not already in RDF
                    # update auto save
                    self.monitor.do.auto_save()
                    self.monitor.do.generate_graph_report()
                    self.monitor.do.send_email_to_users()

                # noinspection PyUnboundLocalVariable
                pbar.update(1)

            IN = 0

        assert IN == 0 # must be out.

        if rAnk == mAster_rank:
            pbar.close()
            print(flush=True)
            # Save it when iterations are done.
            if self.monitor.RDF_filename is not None:

                RDF_filename = self.monitor.RDF_filename

                # noinspection PyBroadException
                try:
                    self.RDF.to_csv(self.monitor.RDF_filename, header=True)

                except:
                    sleep(5) # wait 5 seconds
                    # noinspection PyBroadException
                    try: # try once more
                        self.RDF.to_csv(self.monitor.RDF_filename, header=True)
                    except:
                        sleep(5) # wait 5 seconds
                        # noinspection PyBroadException
                        try: # try once more
                            self.RDF.to_csv(self.monitor.RDF_filename, header=True)
                        except: # save it to a new file!
                            RDF_filename = self.monitor.RDF_filename[:-4] + \
                                           '__completion_save__' + \
                                           randomStringDigits(10) + \
                                           '.csv'

                            self.RDF.to_csv(RDF_filename, header=True)

                # when we save RDF, by default we save the iterator itself only for saving some info (not for resume).
                filename = RDF_filename
                if filename[-4:] == '.csv': filename = filename[:-4]

                if self._save_to_mitr_: # we only do this when we set ``save_to_miter`` to be True.

                    # noinspection PyBroadException
                    try: # try to save the iterator to .mitr file.

                        with open(filename + '.mitr', 'wb') as output:
                            # '.mitr' stands for mimetic iterator.
                            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
                        output.close()

                    except:  # save it to a new file!
                        # noinspection PyBroadException
                        try: # may be the solver property is stopping the saving.
                            self._solver_ = None
                            with open(filename+'.mitr', 'wb') as output:
                                pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
                            output.close()
                        except:
                            pass

    @classmethod
    def read(cls, filename):
        """ """
        cOmm.barrier() # this is important to make the program safe.
        if rAnk == mAster_rank:
            if filename[-4:] == '.csv':
                return pd.read_csv(filename, index_col=0)

            else:

                if filename[-5:] != '.mitr': filename += '.mitr'

                with open(filename, 'rb') as INPUT:
                    iterator = pickle.load(INPUT)
                INPUT.close()
                return iterator
        else:
            # `read` can not resume iterator; just for retrieve info like: t0, dt, steps, solver source, solver dir ...
            return None

    # ==============================================================================================