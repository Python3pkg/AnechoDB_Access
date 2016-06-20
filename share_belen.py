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
### Server
#==============================================================================
class Connection:
    def __init__(self,host):
        ''' 
        Put the host of the database.         
        '''
        self.host=host
        
    def _find_by_var(self,link='',var='',dat = True):
        ''' 
        Return the link or the value of the var entry in the database.
        If dat=True: return link else dat=False return value of var         
        '''
        #dat=True return link, dat=False return id
        l=[]
        r = requests.get(os.path.join(self.host+'/anechodb/api/v1/'+link))
        d = json.loads(r.text)
        s=d['collection']['items']
        
        if dat and var:
            for i in range (np.shape(s)[0]):
                for j in range (np.shape(s[i]['data'])[0]):
                    if var in s[i]['data'][j].values():
                        l=s[i]['href']
                        break
        elif var:
            for i in range (np.shape(s)[0]):
                for j in range (np.shape(s[i]['links'])[0]):
                    if var in s[i]['links'][j]['rel']:
                         l.append(s[i]['data'][0]['value'])
        return l
    def read_link(self,link,idl):
        ''' 
        Print the json collection of the chosen link entry.
        Input:
            link(string):the link to the page of the database. It can be only:
                         'operators','instruments','projects','measurements','beams'.
            idl(int): identifier of the page of the link .        
        '''
        r = requests.get(os.path.join(self.host+'/anechodb/api/v1/'+link+
                        '/%d'%idl))
        if r.status_code==200:
            print(json.loads(r.text))
    
    def search_meas_by_instrument(self,var=''):
        ''' 
        Search wich measurements are linked at the instrument decided by var entry.
        Input:
            var(string):The instrument used for the search (example 'VNA').
        Output:
            m_id(array of int):The identifier of the measurement that use the 
                               instrument.
        '''

        if var:
            m_id=[]
            idl=Connection._find_by_var(self,'instruments',var,True)
            if idl:
                idl=idl.split('/')[-1]
                rel="/api/v1/instruments/"+str(idl)
                m_id=Connection._find_by_var(self,'measurements',rel,False)
            return m_id
        
    def search_meas_by_projects(self,var=''):
        ''' 
        Search wich measurements are linked at the project decided by var entry.
        Input:
            var(string):The project used for the search (example 'LSPE').
        Output:
            m_id(array of int):The identifier of the measurement that use the 
                               project.
        '''

        if var:
            m_id=[]
            idl=Connection._find_by_var(self,'projects',var,True)
            if idl:
                idl=idl.split('/')[-1]
                rel="/api/v1/projects/"+str(idl)
                m_id=Connection._find_by_var(self,'measurements',rel,False)
            return m_id
    def search_beam_by_meas (self,m_id=0):
        ''' 
        Search wich beams are linked at the measurement identifier decided 
        by m_id entry.
        Input:
            m_id(int):The measurement identifier used for the search (example 1).
        Output:
            b_id(array of int):The identifier of the beams linked at the
                               chosen measurement.
        '''        
        
        b_id=[]
        if m_id:
            rel="/api/v1/measurements/"+str(m_id)
            b_id=(Connection._find_by_var(self,'beams',rel,False))
        return b_id
    def get_beam_in_dict_by_id(self,b_id):
        ''' 
        Download the beam chosen by identifier as a dict variable.
        Input:
            b_id(int):The beam identifier to download (example 1).
        Output:
            beam(dict):The beam downloaded. It has 4 fields as the original .h5
                       file plus the attribute field with some extra information.
        '''
        beam={}
        #Connect and dowload the chosen beam
        head = {'Accept': 'application/x-hdf5'}
        r=requests.get(os.path.join(self.host+'/anechodb/api/v1/beams/%d'%b_id),
                       headers=head)
        f_b, p_b = tempfile.mkstemp(suffix='.h5')
        try:
            with os.fdopen(f_b, 'wb') as tmp:
                tmp.write(r.content)
            #create dict variable beam
            fid=h5py.File(p_b,'r')
            for key in fid.keys():
                if key=='DUT' or key=='REF':
                    D=fid['%s'%key]
                    A={}
                    for k in D.keys():
                        P={}
                        P['Amplitude']=D['%s/Ampl'%k].value 
                        P['Phase']=D['%s/Phase'%k].value
                        A['%s'%k]=P    
                    beam['%s'%key]=A
                else:
                    beam['%s'%key]=fid['%s'%key].value
            A={}
            for key in fid.attrs.keys():
                A['%s'%key]=str(fid.attrs.get('%s'%key))
            beam['Attributes']=A   
            fid.close()
        finally:
            os.remove(p_b)
        return beam
#==============================================================================
### Computation
#==============================================================================
class Computation:

    def make_beam_meanvar(beam,f=[],start=0,stop=-1):
        ''' 
        Apply mean and variance at the data stored in beam at the chosen frequencies.
        Input:
            beam(dict):The beam to be computed
            f(array of float):The frequencies of the measure to be computed.
                              If empty all the frequencies in beam will be used.
            start(int):Starting index of the measurement array for the computation.
            stop(int):Stopping index of the measurement array for the computation.
                                            
        Output:
            b(dict):The input beam with measurement changed with the mean and with
                    a new field called Amplitude_Variance with the variance.
        '''
        b=copy.deepcopy(beam)
        if not f:
            f=b['Frequencies']
        for i in range (len(f)):
            id_f=np.where(b['Frequencies']==f[i])
            if b['DUT']['F_%d'%id_f[0]]['Amplitude'].ndim>1:
                b['DUT']['F_%d'%id_f[0]]['Amplitude_Variance']=np.var(b['DUT']
                            ['F_%d'%id_f[0]]['Amplitude'][:,start:stop],axis=1)    
                b['DUT']['F_%d'%id_f[0]]['Amplitude']=np.mean(b['DUT']
                            ['F_%d'%id_f[0]]['Amplitude'][:,start:stop],axis=1)
        return b
        
        
    def center_norm_beam(beam,f=[], center = True, norm = True):
        ''' 
        Apply normalization and centering at the data stored in beam.
        Input:
            beam(dict):The beam to be computed. If Amplitude in beam is a matrix,
                        the mean of the matrix will be used for this computation.
            f(array of float):The frequencies of the measure to be computed.
                              If empty all the frequencies in beam will be used.
            center(bool or int/float):If center=True, apply centering. If it's 
                                      a number, this will be used to correct the 
                                      position.
            norm(bool or int/float):If norm=True, apply normalization. If it's 
                                      a number, this will be used as normalization
                                      factor.
                                            
        Output:
            b(dict):The input beam with Amplitude and Positions computed. The
                    positions of the original beam are stored in Original_Positions
                    field and a new field called Correction is created with 
                    centering and normalization factors stored for each frequency.
        '''
        b=copy.deepcopy(beam)  
        corr={}
        P={}     
        if not f:
            f=b['Frequencies']
        angle=b['Positions'][:,1]    
        for i in range (len(f)):
            c={}
            id_f=np.where(b['Frequencies']==f[i])
    
            if b['DUT']['F_%d'%id_f[0]]['Amplitude'].ndim>1:
                power=np.mean(b['DUT']['F_%d'%id_f[0]]['Amplitude'],axis=1)
            else:
                power=b['DUT']['F_%d'%id_f[0]]['Amplitude']    
                
            
            # Find window at 3 dB
            maxpower = np.max(power)
            index_main_beam = np.where(power >= maxpower -3.)[0]
            main_beam_power = power[index_main_beam]
            main_beam_angle = angle[index_main_beam]
            
            # Interpolate with parabola
            parabola_fit = np.polyfit(main_beam_angle, main_beam_power, 2)
        
            # Find parabola vertex
            vertex_angle = -parabola_fit[1]/2./parabola_fit[0]
            
            # Find parabola maximum
            det = parabola_fit[1]**2 - 4.*parabola_fit[0]*parabola_fit[2]
            vertex_power = -det / (4. * parabola_fit[0])
    
            if type(center == bool):
                if center == True:
                    newangle = angle - vertex_angle
                else:
                    newangle = angle
        
            if type(center) == float or type(center) == int:
                newangle = angle - float(center)            
                    
            if type(norm == bool):
                if norm == True:
                    newpower = power - vertex_power
                else:
                    newpower = power
            
            if type(norm) == int or type(norm) == float:
                newpower = power - norm    
            b['DUT']['F_%d'%id_f[0]]['Amplitude']=newpower
            P['F_%d'%id_f[0]]=newangle
            c['Center']=vertex_angle
            c['Norm']=vertex_power
            corr['F_%d'%id_f[0]]=c
        b['Original_positions']=b['Positions']
        b['Positions']=P
        b['Correction']=corr
        return b
        
    def phase_center(L,d,angle,magn):
        # d:distance between receiver feed(phase center) and rotation centre[m]
        # L:distance between receiver feed and trasmitter feed[m] 
    
        corr_ang=np.rad2deg(np.arctan2(d*np.sin(np.deg2rad(angle)),L+d*(1-np.cos(np.deg2rad(angle)))))
        newangle=angle+corr_ang
        L_true=np.sqrt((d*np.sin(np.deg2rad(angle)))**2+(L+d*(1-np.cos(np.deg2rad(angle))))**2)
        corr_pow=20*np.log(L_true/(L+d))
        newpower=magn+corr_pow
        corr_phase=(corr_ang,corr_pow)
        return (corr_phase,newangle, newpower)
        
    def sim_diff (angle_sim,angle_mis,magn_sim,magn_mis):
        tot_ind=np.zeros(np.shape(angle_mis))
        diff=np.zeros(np.shape(angle_mis))
        for i in range(0, np.shape(angle_mis)[0]):
            index=np.where(angle_sim>=angle_mis[i])[0][0]
            if ((angle_mis[i]-angle_sim[index-1]) < (angle_sim[index]-angle_mis[i])):
                index=index-1
            else:
                pass
            tot_ind[i]=index
            diff[i]=magn_mis[i]-magn_sim[index]
            del index
        return (diff,tot_ind)