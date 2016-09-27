# -*- coding: utf-8 -*-
"""
@author: Jacopo Martelli
"""
import numpy as np
import copy
'''
Module with some useful function to compute simple analytics correction to beam patterns.
'''


def make_beam_meanvar(
        beam: dict,
        f: list=[],
        start: int=0,
        stop: int=None)->dict:
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
        if (b['DUT']['F_%d' % id_f[0][0]]['Ampl'].ndim > 1) and id_f:
            b['DUT'][
                'F_%d' %
                id_f[0][0]]['Ampl_Var'] = np.var(
                b['DUT'][
                    'F_%d' %
                    id_f[0][0]]['Ampl'][
                    :,
                    start:stop],
                axis=1)
            b['DUT'][
                'F_%d' %
                id_f[0][0]]['Ampl'] = np.mean(
                b['DUT'][
                    'F_%d' %
                    id_f[0][0]]['Ampl'][
                    :,
                    start:stop],
                axis=1)
    return b


def center_norm_beam(beam: dict, f: list=[], center=True, norm=True)->dict:
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
        raise Exception(
            'Input beam is a crosspolar, so center and norm entries must be float or int')
    else:
        if not isinstance(f, list):
            f = list(f)
        if not f:
            f = b['Frequencies']
        angle = b['Positions'][:, 1]
        for i in range(len(f)):
            c = {}
            id_f = np.where(b['Frequencies'] == f[i])

            if b['DUT']['F_%d' % id_f[0][0]]['Ampl'].ndim > 1:
                power = np.mean(
                    b['DUT'][
                        'F_%d' %
                        id_f[0][0]]['Ampl'],
                    axis=1)
            else:
                power = b['DUT']['F_%d' % id_f[0][0]]['Ampl']

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

            if type(center == bool):
                if center:
                    newangle = angle - vertex_angle
                else:
                    newangle = angle

            if type(center) == float or type(center) == int:
                vertex_angle = float(center)
                newangle = angle - vertex_angle

            if type(norm == bool):
                if norm:
                    newpower = power - vertex_power
                else:
                    newpower = power

            if type(norm) == int or type(norm) == float:
                vertex_power = float(norm)
                newpower = power - vertex_power

            b['DUT']['F_%d' % id_f[0][0]]['Ampl'] = newpower
            P['F_%d' % id_f[0][0]] = newangle
            c['Center'] = vertex_angle
            c['Norm'] = vertex_power
            corr['F_%d' % id_f[0][0]] = c
        b['Original_positions'] = b['Positions']
        b['Positions'] = P
        b['Correction'] = corr
    return b
