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
        self.host=host
        
    def find_by_var(self,link='',var='',dat = True):
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
    def read_link(self,link='',idl=''):
        r = requests.get(os.path.join(self.host+'/anechodb/api/v1/'+link+'/'+idl))
        if r.status_code==200:
            print(json.loads(r.text))
    
    def search_meas_by_instrument(var=''):
        if var:
            m_id=[]
            idl=Connection.find_by_var('instruments',var,True)
            if idl:
                idl=idl.split('/')[-1]
                rel="/api/v1/instruments/"+str(idl)
                m_id=Connection.find_by_var('measurements',rel,False)
            return m_id
        
    def search_meas_by_projects(var=''):
        if var:
            m_id=[]
            idl=Connection.find_by_var('projects',var,True)
            if idl:
                idl=idl.split('/')[-1]
                rel="/api/v1/projects/"+str(idl)
                m_id=Connection.find_by_var('measurements',rel,False)
            return m_id
    def search_beam_by_meas (m_id=0):
        b_id=[]
        if m_id:
            rel="/api/v1/measurements/"+str(m_id)
            b_id=(Connection.find_by_var('beams',rel,False))
        return b_id
    def get_beam_in_dict_by_id(self,b_id):
        beam={}
        #Connect and dowload the chosen beam
        head = {'Accept': 'application/x-hdf5'}
        r=requests.get(os.path.join(self.host+'/anechodb/api/v1/beams/%s'%b_id),headers=head)
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
class computation:

    def make_beam_meanvar(beam,f=[],start=0,stop=-1):
        b=copy.deepcopy(beam)
        if not f:
            f=b['Frequencies']
        for i in range (len(f)):
            id_f=np.where(b['Frequencies']==f[i])
            if b['DUT']['F_%d'%id_f[0]]['Amplitude'].ndim>1:
                b['DUT']['F_%d'%id_f[0]]['Amplitude_Variance']=np.var(b['DUT']['F_%d'%id_f[0]]['Amplitude'][:,start:stop],axis=1)    
                b['DUT']['F_%d'%id_f[0]]['Amplitude']=np.mean(b['DUT']['F_%d'%id_f[0]]['Amplitude'][:,start:stop],axis=1)
        return b
        
        
    def center_norm_beam(beam,f=[], center = True, norm = True):
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