# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 11:44:55 2016

@author: Jacopo Martelli
"""

import numpy as np
import requests
import json
import os
import tempfile
import h5py


class Connection:
    '''
    Class used for the comunication with the database to find and retrieve beams stored in it.
    '''

    def __init__(self, host):
        '''
        Put the host of the database.
        '''
        self.host = host

    def _find_by_var(self, link: str='', var: str='', dat: bool=True)->list:
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
                        break  # there should be only a link for each entry var
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

    def search_meas_by_instruments(self, var: str='')->list:
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
                raise Exception('Nothing found with this name: %s' % var)

    def search_meas_by_projects(self, var: str='')->list:
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
                raise Exception('Nothing found with this name: %s' % var)

    def search_beam_by_meas(self, m_id: int=0)->list:
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
            raise Exception('No beam linked to measurement id: %d' % m_id)

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
