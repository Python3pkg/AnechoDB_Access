# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 10:52:46 2016

@author: Jacopo Martelli
"""
import numpy as np
import requests
import copy
import json
import os
import tempfile
import h5py
#==============================================================================
# Server
#==============================================================================


class Connection:
    '''
    Class used for the comunication with the database to find and retrieve beams stored in it.
    '''

    def __init__(self, host):
        '''
        Put the host of the database.
        '''
        self.host = host

    def _find_by_var(self, link:str='', var:str='', dat:bool=True)->list:
        '''
        Return the link or the value of the var entry in the database.

        If dat=True return link else, dat=False, return value of var.
        '''

        l = []
        r = requests.get(os.path.join(self.host + '/anechodb/api/v1/' + link))
        d = json.loads(r.text)
        s = d['collection']['items']

        if dat and var:
            for i in range(np.shape(s)[0]):
                for j in range(np.shape(s[i]['data'])[0]):
                    if var in s[i]['data'][j].values():
                        l = s[i]['href']
                        break#there should be only a link for each entry var
        elif var:
            for i in range(np.shape(s)[0]):
                for j in range(np.shape(s[i]['links'])[0]):
                    if var in s[i]['links'][j]['rel']:
                        l.append(s[i]['data'][0]['value'])
        return l

    def print_link(self, link, idl):
        '''
        Print the json collection of the chosen link entry.

        Input:

            link(string):the link to the page of the database. It can only be:
                         'operators','instruments','projects','measurements','beams'.

            idl(int): identifier of the page of the link .
        '''
        r = requests.get(os.path.join(self.host + '/anechodb/api/v1/' + link +
                                      '/%d' % idl))
        if r.status_code == 200:
            print(json.loads(r.text))

    def search_meas_by_instruments(self, var:str='')->list:
        '''
        Search which measurements are linked at the instrument decided by var entry.

        Input:

            var(string):The instrument used for the search (example 'VNA').

        Output:

            m_id(list of int):The identifier of the measurement that use the instrument.
        '''

        if var:
            m_id = []
            idl = Connection._find_by_var(self, 'instruments', var, True)
            if idl:
                idl = idl.split('/')[-1]
                rel = "/api/v1/instruments/" + str(idl)
                m_id = Connection._find_by_var(
                    self, 'measurements', rel, False)
                return m_id
            else:
                raise Exception('Nothing found with this name: %s'%var)

    def search_meas_by_projects(self, var:str='')->list:
        '''
        Search which measurements are linked at the project decided by var entry.

        Input:

            var(string):The project used for the search (example 'LSPE').

        Output:

            m_id(list of int):The identifier of the measurement that use the project.
        '''

        if var:
            m_id = []
            idl = Connection._find_by_var(self, 'projects', var, True)
            if idl:
                idl = idl.split('/')[-1]
                rel = "/api/v1/projects/" + str(idl)
                m_id = Connection._find_by_var(
                    self, 'measurements', rel, False)
                return m_id
            else:
                raise Exception('Nothing found with this name: %s'%var)

    def search_beam_by_meas(self, m_id:int=0)->list:
        '''
        Search which beams are linked at the measurement identifier decided by m_id entry.

        Input:

            m_id(int):The measurement identifier used for the search (example 1).

        Output:

            b_id(list of int):The identifier of the beams linked at the chosen measurement.
        '''

        b_id = []
        if m_id:
            rel = "/api/v1/measurements/" + str(m_id)
            b_id = (Connection._find_by_var(self, 'beams', rel, False))
            return b_id
        else:
            raise Exception('No beam linked to measurement id: %d'%m_id)

    def get_beam_in_dict_by_id(self, b_id: int) -> dict:
        '''
        Download the beam chosen by identifier as a dict variable.

        Input:

            b_id(int):The beam identifier to download (example 1).

        Output:

            beam(dict):The beam downloaded. It has 4 fields as the original .h5 file plus the attribute field with some extra information.
        '''
        beam = {}
        # Connect and dowload the chosen beam
        head = {'Accept': 'application/x-hdf5'}
        r = requests.get(
            os.path.join(
                self.host +
                '/anechodb/api/v1/beams/%d' %
                b_id),
            headers=head)
        f_b, p_b = tempfile.mkstemp(suffix='.h5')
        try:
            with os.fdopen(f_b, 'wb') as tmp:
                tmp.write(r.content)
            # Create dict variable beam
            fid = h5py.File(p_b, 'r')
            for key in fid.keys():
                if key == 'DUT' or key == 'REF':
                    D = fid['%s' % key]
                    A = {}
                    for k in D.keys():
                        P = {}
                        P['Amplitude'] = D['%s/Ampl' % k].value
                        P['Phase'] = D['%s/Phase' % k].value
                        A['%s' % k] = P
                    beam['%s' % key] = A
                else:
                    beam['%s' % key] = fid['%s' % key].value
            A = {}
            for key in fid.attrs.keys():
                A['%s' % key] = str(fid.attrs.get('%s' % key))
            beam['Attributes'] = A
            fid.close()
        finally:
            os.remove(p_b)
        return beam
#==============================================================================
# Computation
#==============================================================================


class Computation:
    '''
    Class with various function useful to apply corrections at the beam patterns.
    '''
    def make_beam_meanvar(beam:dict, f:list=[], start:int=0, stop:int=None)->dict:
        '''
        Apply mean and variance at the data stored in beam at the chosen frequencies.

        Input:

            beam(dict):The beam to be computed. The Amplitude field stored should be a matrix to apply mean and variance along the measurement points.

            f(list of float):The frequencies of the measure to be computed. If empty all the frequencies in beam are used. If it's only a number use as input a list(Example f=[40])

            start(int):Starting index of the measurement array for the computation.

            stop(int):Stopping index of the measurement array for the computation.

        Output:

            b(dict):The input beam with measurement changed with the mean and with a new field called Amplitude_Variance with the variance.
        '''
        b = copy.deepcopy(beam)
        if not isinstance(f, list):
            f = list(f)
        if not f:
            f = b['Frequencies']
        for i in range(len(f)):
            id_f = np.where(b['Frequencies'] == f[i])
            if (b['DUT']['F_%d' % id_f[0][0]]['Amplitude'].ndim > 1) and id_f:
                b['DUT'][
                    'F_%d' %
                    id_f[0][0]]['Amplitude_Variance'] = np.var(
                    b['DUT'][
                        'F_%d' %
                        id_f[0][0]]['Amplitude'][
                        :,
                        start:stop],
                    axis=1)
                b['DUT'][
                    'F_%d' %
                    id_f[0][0]]['Amplitude'] = np.mean(
                    b['DUT'][
                        'F_%d' %
                        id_f[0][0]]['Amplitude'][
                        :,
                        start:stop],
                    axis=1)
        return b

    def center_norm_beam(beam:dict, f:list=[], center=True, norm=True)->dict:
        '''
        Apply normalization and centering at the data stored in beam.

        Input:

            beam(dict):The beam to be computed. If Amplitude in beam is a matrix, the mean of the matrix is used for this computation.

            f(array of float):The frequencies of the measure to be computed. If empty all the frequencies in beam are used. If it's only a number use as input a list(Example f=[40])

            center(bool or int/float):If center=True, apply centering. If it's a number, this is used to correct the position.

            norm(bool or int/float):If norm=True, apply normalization. If it's a number, this is used as normalization factor.

        Output:

            b(dict):The input beam with Amplitude and Positions computed. The positions of the original beam are stored in Original_Positions
            field and a new field called Correction is created with centering and normalization factors stored for each frequency.

        Notes:

            If the beam is not copolar (it's seen in the Attributes field) input variables center and norm MUST be numbers to use this function.
        '''
        b = copy.deepcopy(beam)
        corr = {}
        P = {}
        if all([b['Attributes']['Type'][-1] != 'O',
                (isinstance(center, bool) or isinstance(norm, bool))]):
            raise Exception('Input beam is a crosspolar, so center and norm entries must be float or int')
        else:
            if not isinstance(f, list):
                f = list(f)
            if not f:
                f = b['Frequencies']
            angle = b['Positions'][:, 1]
            for i in range(len(f)):
                c = {}
                id_f = np.where(b['Frequencies'] == f[i])

                if b['DUT']['F_%d' % id_f[0][0]]['Amplitude'].ndim > 1:
                    power = np.mean(
                        b['DUT'][
                            'F_%d' %
                            id_f[0][0]]['Amplitude'],
                        axis=1)
                else:
                    power = b['DUT']['F_%d' % id_f[0][0]]['Amplitude']

                # Find window at 3 dB
                maxpower = np.max(power)
                index_main_beam = np.where(power >= maxpower - 3.)[0]
                main_beam_power = power[index_main_beam]
                main_beam_angle = angle[index_main_beam]

                # Interpolate with parabola
                parabola_fit = np.polyfit(main_beam_angle, main_beam_power, 2)

                # Find parabola vertex
                vertex_angle = -parabola_fit[1] / 2. / parabola_fit[0]

                # Find parabola maximum
                det = parabola_fit[1]**2 - 4. * \
                    parabola_fit[0] * parabola_fit[2]
                vertex_power = -det / (4. * parabola_fit[0])

                if isinstance(center, bool):
                    if center:
                        newangle = angle - vertex_angle
                    else:
                        newangle = angle

                if isinstance(center, float) or isinstance(center, int):
                    vertex_angle = float(center)
                    newangle = angle - vertex_angle

                if isinstance(norm, bool):
                    if norm:
                        newpower = power - vertex_power
                    else:
                        newpower = power

                if isinstance(norm, int) or isinstance(norm, float):
                    vertex_power = float(norm)
                    newpower = power - vertex_power

                b['DUT']['F_%d' % id_f[0][0]]['Amplitude'] = newpower
                P['F_%d' % id_f[0][0]] = newangle
                c['Center'] = vertex_angle
                c['Norm'] = vertex_power
                corr['F_%d' % id_f[0][0]] = c
            b['Original_positions'] = b['Positions']
            b['Positions'] = P
            b['Correction'] = corr
        return b

    def phase_center(L:float, d:float, angle, magn):
        '''
         Compute correction in the case the phase center of the antenna is not placed at the rotation axis of the measuring system.

        Input:

            L(int or float):distance between receiver feed and transmitter feed[m].

            d(int or float):distance between phase center of receiver feed and rotation centre[m].

            angle(array):The positions of the measure.

            magn(array):The amplitude of the measure.

        Output:

            corr_phase(tuple of float): Containing the correction factors applied at the positions and at the magnitudes.

            newangle(array): The positions with the correction applied.

            newpower(array): The amplitudes with the correction applied.
        '''

        corr_ang = np.rad2deg(np.arctan2(
            d * np.sin(np.deg2rad(angle)), L + d * (1 - np.cos(np.deg2rad(angle)))))
        newangle = angle + corr_ang
        L_true = np.sqrt((d * np.sin(np.deg2rad(angle)))**2 +
                         (L + d * (1 - np.cos(np.deg2rad(angle))))**2)
        corr_pow = 20 * np.log(L_true / (L + d))
        newpower = magn + corr_pow
        corr_phase = (corr_ang, corr_pow)
        return (corr_phase, newangle, newpower)

    def sim_diff(angle_sim, angle_mis, magn_sim, magn_mis):
        '''
        Compute difference between simulation and measured data.
        Check the nearest angle between simulation and measure and apply difference between them.
        More discrete are the simulation angles more the computation is correct.

        N.B. It's up to the user to choose the same frequency for simulation and measure.

        Input:

            angle_sim(array):The positions used in simulation.

            angle_mis(array):The positions of the measure.

            magn_sim(array):The amplitude obtained from simulation.

            magn_mis(array):The amplitude of the measure.

        Output:

            diff(array):The array containing the difference computed.
        '''
        diff = np.zeros(np.shape(angle_mis))
        for i in range(0, np.shape(angle_mis)[0]):
            index = np.where(angle_sim >= angle_mis[i])[0][0]
            if ((angle_mis[i] - angle_sim[index - 1])
                    < (angle_sim[index] - angle_mis[i])):
                index = index - 1
            else:
                pass
            diff[i] = magn_mis[i] - magn_sim[index]
            del index
        return diff
