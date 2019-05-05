"""
Module for setting up some commonly used fields for MAESTRO urca datasets with yt.

Donald E. Willcox
"""

import yt
import numpy as np

class PhysicalConstants:
    N_AVO = 6.02214129e23
    
class DatasetHelpers:
    @staticmethod
    def get_field(ds, field_name):
        field = None
        field_short_name = None
        for f in ds.field_list + ds.derived_field_list:
            if f[1] == field_name:
                field_short_name = f[1]
                field = f
                return field, field_short_name
        if not field:
            print('Field {} not present.'.format(field_name))
            return None, None

class UrcaShellFields(object):
    def __init__(self):
        return

    def setup(self, ds):
        # ds should be a MAESTRO dataset in yt corresponding to the urca shell variables

        try:
            ds.add_field(('boxlib','urca23_shell_unscaled'), sampling_type='cell', units='',
                         function=UrcaShellFields._urca23_shell_unscaled)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: urca23_shell_unscaled field could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','urca23_shell'), sampling_type='cell', units='',
                         function=UrcaShellFields._urca23_shell)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: urca23_shell field could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','weak_xrate_na23'), sampling_type='cell', units='',
                         function=UrcaShellFields._weak_xrate_na23)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: weak_xrate_na23 field could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','weak_xrate_ne23'), sampling_type='cell', units='',
                         function=UrcaShellFields._weak_xrate_ne23)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: weak_xrate_ne23 field could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','reduced_x23'), sampling_type='cell', units='',
                         function=UrcaShellFields._reduced_x23)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: reduced_x23 could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','enucloss_epart_urca23'), sampling_type='cell', units='g',
                         function=UrcaShellFields._enucloss_epart_urca23)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: enucloss_epart_urca23 could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','enucdot_dqweak_urca23'), sampling_type='cell', units='g',
                         function=UrcaShellFields._enucdot_dqweak_urca23)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: enucdot_dqweak_urca23 could not be added because it relies on a field not in the dataset.')
            pass
        
        try:
            ds.add_field(('boxlib','enucloss_sneut'), sampling_type='cell', units='g',
                         function=UrcaShellFields._enucloss_sneut)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: enucloss_sneut could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','enucdot_ion_binding'), sampling_type='cell', units='g',
                         function=UrcaShellFields._enucdot_ion_binding)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: enucdot_ion_binding could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','enucdot_summed_total'), sampling_type='cell', units='g',
                         function=UrcaShellFields._enucdot_summed_total)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: enucdot_summed_total could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','enucdot_total'), sampling_type='cell', units='erg/s',
                         function=UrcaShellFields._enucdot_total)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: enucdot_total could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','sum_omegadots'), sampling_type='cell', units='1/s',
                         function=UrcaShellFields._sum_omegadots)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: sum_omegadots could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','sum_omegadot_urca23'), sampling_type='cell', units='1/s',
                         function=UrcaShellFields._sum_omegadot_urca23)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: sum_omegadot_urca23 could not be added because it relies on a field not in the dataset.')
            pass

        try:
            ds.add_field(('boxlib','xc12_complement'), sampling_type='cell', units='',
                         function=UrcaShellFields._xc12_complement)
        except yt.utilities.exceptions.YTFieldNotFound:
            print('WARNING: xc12_complement could not be added because it relies on a field not in the dataset.')
            pass
        
    @staticmethod
    def _urca23_shell_unscaled(field, data):
        return data['boxlib','ecap23']*data['boxlib','beta23']*data['boxlib','X(na23)']*data['boxlib','X(ne23)']

    @staticmethod
    def _urca23_shell(field, data):
        return data['boxlib','urca23_shell_unscaled']/np.amax(data['boxlib','urca23_shell_unscaled'])

    @staticmethod
    def _weak_xrate_na23(field, data):
        return data['boxlib','beta23']*data['boxlib','X(ne23)'] - data['boxlib', 'X(na23)']*data['boxlib','ecap23']

    @staticmethod
    def _weak_xrate_ne23(field, data):
        return data['boxlib', 'X(na23)']*data['boxlib','ecap23'] - data['boxlib','beta23']*data['boxlib','X(ne23)']

    @staticmethod
    def _reduced_x23(field, data):
        return data['boxlib','X(na23)']*data['boxlib', 'X(ne23)']/(data['boxlib','X(na23)']+data['boxlib', 'X(ne23)'])

    @staticmethod
    def _enucloss_epart_urca23(field, data):
        return (-data['boxlib', 'epart_ecap23'] - data['boxlib', 'epart_beta23']) * data['boxlib', 'density'] * data['boxlib', 'cell_volume']

    @staticmethod
    def _enucloss_sneut(field, data):
        # Energy loss rate due to plasma and other neutrino losses (not Urca) in erg/g/s * g/cm^3 * cm^3 = erg/s
        return data['boxlib', 'sneut'] * data['boxlib', 'density'] * data['boxlib', 'cell_volume']

    @staticmethod
    def _enucdot_ion_binding(field, data):
        # Energy generation rate due to ion binding energies (does not include dQ corrections for Urca reactions) in erg/g/s * g/cm^3 * cm^3 = erg/s
        return data['boxlib', 'ionenuc'] * data['boxlib', 'density'] * data['boxlib', 'cell_volume']

    @staticmethod
    def _enucdot_dqweak_urca23(field, data):
        return (data['boxlib', 'dqweak_ecap23'] + data['boxlib', 'dqweak_beta23']) * data['boxlib', 'density'] * data['boxlib', 'cell_volume']

    @staticmethod
    def _enucdot_summed_total(field, data):
        return (data['boxlib', 'enucdot_ion_binding'] - data['boxlib', 'enucloss_sneut'] - data['boxlib', 'enucloss_epart_urca23'] + data['boxlib', 'enucdot_dqweak_urca23'])

    @staticmethod
    def _enucdot_total(field, data):
        return data['boxlib', 'enucdot'] * data['boxlib', 'density'] * data['boxlib', 'cell_volume'].in_units('cm**3')

    @staticmethod
    def _sum_omegadots(field, data):
        return np.sum(data['microphysics', 'omegadots'], axis=0)

    @staticmethod
    def _sum_omegadot_urca23(field, data):
        return data['boxlib', 'omegadot(ne23)'] + data['boxlib', 'omegadot(na23)']

    @staticmethod
    def _xc12_complement(field, data):
        return 0.5 - data['boxlib', 'X(c12)']
