# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 00:12:35 2018

@author: HP
"""

import numpy as np
import matplotlib.pyplot as pl

phase = 'Ca--HG'
pressure = '2410'
serial = [100,300,600,900,1200,2000,3000] # temperature range

filename = 'XDATCAR.' + str(serial[0])
XDATCAR = open(filename,'r')
parax = []
paray = []
paraz = []
ratio = 0
system_name = ''
atom_name = ''
atom_num = 0
frame = 0
count = 0
for line in XDATCAR:
    count += 1
    line.strip('\n')
    line = line.split()
    if count == 1:
        system_name = line[0]
    elif count == 2:
        ratio = float(line[0])
    elif count == 3:
        parax.append(float(line[0]))
        parax.append(float(line[1]))
        parax.append(float(line[2]))
    elif count == 4:
        paray.append(float(line[0]))
        paray.append(float(line[1]))
        paray.append(float(line[2]))
    elif count == 5:
        paraz.append(float(line[0]))
        paraz.append(float(line[1]))
        paraz.append(float(line[2]))
    elif count == 6:
        atom_name = line[0]
    elif count == 7:
        atom_num = int(line[0])
    elif line[0] == 'Direct':
        frame += 1
XDATCAR.close()
    
frame = 10000
sta_num = 3001
sta_time = 1001

para_x = parax[0]
para_y = paray[1]
para_z = paraz[2]
posx = np.zeros([sta_num,atom_num])
posy = np.zeros([sta_num,atom_num])
posz = np.zeros([sta_num,atom_num])

def MSD(posx,posy,posz,steps):
    POSX = []
    POSY = []
    POSZ = []
    posx_ = 0
    posy_ = 0
    posz_ = 0
    for i in range(steps - 1):
        posx__ = posx[i + 1] - posx[i]
        posx__ = (posx__ > 0.5) * (-1.0 + posx__) + \
            (-posx__ > 0.5) * (1.0 + posx__) + (abs(posx__) < 0.5) * posx__
        posy__ = posy[i + 1] - posy[i]
        posy__ = (posy__ > 0.5) * (-1.0 + posy__) + \
            (-posy__ > 0.5) * (1.0 + posy__) + (abs(posy__) < 0.5) * posy__
        posz__ = posz[i + 1] - posz[i]
        posz__ = (posz__ > 0.5) * (-1.0 + posz__) + \
            (-posz__ > 0.5) * (1.0 + posz__) + (abs(posz__) < 0.5) * posz__
        posx_ += posx__
        posy_ += posy__
        posz_ += posz__
        POSX.append((posx_ * para_x) ** 2)
        POSY.append((posy_ * para_y) ** 2)
        POSZ.append((posz_ * para_z) ** 2)
    return POSX,POSY,POSZ

def Auto_correlation(posx,posy,posz,steps,time):
    POSX = []
    POSY = []
    POSZ = []
    for i in range(steps - 1):
        posx__ = posx[i + 1] - posx[i]
        posx__ = (posx__ > 0.5) * (-1.0 + posx__) + \
            (-posx__ > 0.5) * (1.0 + posx__) + (abs(posx__) < 0.5) * posx__
        posy__ = posy[i + 1] - posy[i]
        posy__ = (posy__ > 0.5) * (-1.0 + posy__) + \
            (-posy__ > 0.5) * (1.0 + posy__) + (abs(posy__) < 0.5) * posy__
        posz__ = posz[i + 1] - posz[i]
        posz__ = (posz__ > 0.5) * (-1.0 + posz__) + \
            (-posz__ > 0.5) * (1.0 + posz__) + (abs(posz__) < 0.5) * posz__
        POSX.append(posx__ * para_x)
        POSY.append(posy__ * para_y)
        POSZ.append(posz__ * para_z)
    VELX = []
    VELY = []
    VELZ = []
    for i in range(time):
        velx = 0
        vely = 0
        velz = 0
        for j in range(steps - time):
            velx += POSX[j] * POSX[j + i] / (steps - time)
            vely += POSY[j] * POSY[j + i] / (steps - time)
            velz += POSZ[j] * POSZ[j + i] / (steps - time)
        VELX.append(velx)
        VELY.append(vely)
        VELZ.append(velz)
    return VELX,VELY,VELZ

guest_number = 32
guest_list = [18,45,23,43,28,40,33,30,27,39,24,44,17,46,29,34,\
              42,21,37,26,19,47,31,35,32,36,48,20,41,25,38,22]
guest_list = (np.array(guest_list) - 1).tolist()

for k in range(len(serial)):
    line_count = 0
    fra_count = 0
    count_sta = 0
    count_atom = 0
    filename = 'XDATCAR.' + str(serial[k])
    XDATCAR = open(filename,'r')
    for line in XDATCAR:
        line_count += 1
        line.strip('\n')
        line = line.split()
        if line_count >= 8 and line[0] == 'Direct':
            fra_count += 1
            count_atom = 0
        elif fra_count > frame:
            break
        elif fra_count >= frame - sta_num + 1 and line[0] != 'Direct':
            posx[fra_count - frame + sta_num - 1,count_atom] = float(line[0])
            posy[fra_count - frame + sta_num - 1,count_atom] = float(line[1])
            posz[fra_count - frame + sta_num - 1,count_atom] = float(line[2])
            count_atom += 1
    XDATCAR.close()
    print('XDATCAR data inputed...')
    pl.figure('MSD temperature = %d'%serial[k],figsize = [12,8])
    POSX,POSY,POSZ = MSD(posx,posy,posz,sta_num)
    print('MSD calculation done...')
    GUEST_POSX,GUEST_POSY,GUEST_POSZ = [],[],[]
    HOST_POSX,HOST_POSY,HOST_POSZ = [],[],[]
    for i in range(sta_num - 1):
        GUEST_POSX.append(np.mean(POSX[i][guest_list]))
        GUEST_POSY.append(np.mean(POSY[i][guest_list]))
        GUEST_POSZ.append(np.mean(POSZ[i][guest_list]))
        HOST_POSX.append((np.sum(POSX[i]) - np.sum(POSX[i][guest_list]))\
                          /(atom_num - guest_number))
        HOST_POSY.append((np.sum(POSY[i]) - np.sum(POSY[i][guest_list]))\
                          /(atom_num - guest_number))
        HOST_POSZ.append((np.sum(POSZ[i]) - np.sum(POSZ[i][guest_list]))\
                          /(atom_num - guest_number))
    pl.plot(np.arange(1,sta_num),GUEST_POSX,lw = 3,label = 'GUEST_MSD_X')
    pl.plot(np.arange(1,sta_num),GUEST_POSY,lw = 3,label = 'GUEST_MSD_Y')
    pl.plot(np.arange(1,sta_num),GUEST_POSZ,lw = 3,label = 'GUEST_MSD_Z')
    pl.plot(np.arange(1,sta_num),HOST_POSX,lw = 3,label = 'HOST_MSD_X')
    pl.plot(np.arange(1,sta_num),HOST_POSY,lw = 3,label = 'HOST_MSD_Y')
    pl.plot(np.arange(1,sta_num),HOST_POSZ,lw = 3,label = 'HOST_MSD_Z')
    #pl.plot(np.arange(1,sta_num),GUEST_POSZ,lw = 3,label = '%d K'%serial[k])
    pl.title('Ca--HG  %d GPa  %d K'%(int(int(pressure) / 10),serial[k]),fontsize = 18)
    pl.xticks(fontsize = 16)
    pl.xlabel('Time ps',fontsize = 16)
    pl.yticks(fontsize = 16)
    pl.ylabel('MSD A^2',fontsize = 16)
    pl.legend(fontsize = 16,loc = 'best')
    VELX,VELY,VELZ = Auto_correlation(posx,posy,posz,sta_num,sta_time)
    print('Velocity auto-correlation calculation done...')
    pl.figure('VAC temperature = %d'%serial[k],figsize = [12,8])
    GUEST_VELX,GUEST_VELY,GUEST_VELZ = [],[],[]
    HOST_VELX,HOST_VELY,HOST_VELZ = [],[],[]
    for i in range(sta_time):
        GUEST_VELX.append(np.mean(VELX[i][guest_list]))
        GUEST_VELY.append(np.mean(VELY[i][guest_list]))
        GUEST_VELZ.append(np.mean(VELZ[i][guest_list]))
        HOST_VELX.append((np.sum(VELX[i]) - np.sum(VELX[i][guest_list]))\
                          /(atom_num - guest_number))
        HOST_VELY.append((np.sum(VELY[i]) - np.sum(VELY[i][guest_list]))\
                          /(atom_num - guest_number))
        HOST_VELZ.append((np.sum(VELZ[i]) - np.sum(VELZ[i][guest_list]))\
                          /(atom_num - guest_number))
    frequency_range = np.linspace(-500,500,sta_time)
    window = np.arange(int((sta_time - 1) / 2),int((sta_time - 1) / 2) + 51)
    pl.plot(frequency_range[window],abs(np.fft.fftshift(np.fft.fft(GUEST_VELX)))[window]**2,\
            lw = 3,label = 'GUEST_VAC_X')
    pl.plot(frequency_range[window],abs(np.fft.fftshift(np.fft.fft(GUEST_VELY)))[window]**2,\
            lw = 3,label = 'GUEST_VAC_Y')
    pl.plot(frequency_range[window],abs(np.fft.fftshift(np.fft.fft(GUEST_VELZ)))[window]**2,\
            lw = 3,label = 'GUEST_VAC_Z')
    pl.plot(frequency_range[window],abs(np.fft.fftshift(np.fft.fft(HOST_VELX)))[window]**2,\
            lw = 3,label = 'HOST_VAC_X')
    pl.plot(frequency_range[window],abs(np.fft.fftshift(np.fft.fft(HOST_VELY)))[window]**2,\
            lw = 3,label = 'HOST_VAC_Y')
    pl.plot(frequency_range[window],abs(np.fft.fftshift(np.fft.fft(HOST_VELZ)))[window]**2,\
            lw = 3,label = 'HOST_VAC_Z')
    pl.title('Ca--HG  %d GPa  %d K'%(int(int(pressure) / 10),serial[k]),fontsize = 18)
    #pl.plot(frequency_range[window],abs(np.fft.fftshift(np.fft.fft(GUEST_VELZ)))[window]**2,\
    #        lw = 3,label = '%d K'%serial[k])
    pl.title('Ca--HG  %d GPa  GUEST_Z'%int(int(pressure) / 10),fontsize = 18)
    pl.xticks(fontsize = 16)
    pl.xlabel('Frequency THz',fontsize = 16)
    pl.yticks(fontsize = 16)
    pl.ylabel('Power',fontsize = 16)
    pl.legend(fontsize = 16,loc = 'best')