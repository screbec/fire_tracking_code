# -*- coding: utf-8 -*-
""" FireTime
This module include functions used to handle the time conversion
The time step in the tracking system is defined as a list (year, month, day, ampm).
The following functions are used to convert times between different formats.
    t : time steps, tuple (year,month,day,ampm)
    d : date, datetime.date()
    ampm : ampm, str()
    dt : time steps, datetime.datetime()
"""

def t_nb(t,nb='next'):
    ''' Calculate the next or previous time step (year, month, day, ampm)
    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' for present time
    nb : str, 'next'|'previous'
        option to extract next or previous time step

    Returns
    -------
    t_out : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' for next/previous time
    '''
    from datetime import date, timedelta

    # the next time step
    if nb == 'next':
        # if current time is 'AM', set next time as the current day and 'PM'
        if t[-1] == 'AM':
            t_out = list(t[:-1])
            t_out.append('PM')
        # if current time is 'PM', set next time as the following day and 'AM'
        else:
            d = date(*t[:-1])
            d_out = d + timedelta(days=1)
            t_out = [d_out.year,d_out.month,d_out.day,'AM']

    # the previous time step
    elif nb == 'previous':
        # if current time is 'PM', set previous time as the current day and 'AM'
        if t[-1] == 'PM':
            t_out = list(t[:-1])
            t_out.append('AM')
        # if current time is 'AM', set previous time as the previous day and 'PM'
        else:
            d = date(*t[:-1])
            d_out = d + timedelta(days=-1)
            t_out = [d_out.year,d_out.month,d_out.day,'PM']
    return t_out

def t_dif(t1,t2):
    ''' calculate the time difference between two time steps
    Parameters
    ----------
    t1 : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' for time 1
    t2 : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' for time 2

    Returns
    -------
    dt : float
        time difference in days (t2-t1), half day as 0.5
    '''
    from datetime import date

    # calculate the day difference
    d1 = date(*t1[:-1])
    d2 = date(*t2[:-1])
    dt = (d2-d1).days

    # adjust according to ampm difference
    if t1[-1] != t2[-1]:
        if t1[-1] == 'PM':
            dt -= 0.5
        else:
            dt += 0.5
    return dt

def t2d(t):
    ''' convert a t tuple to date and ampm
    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'

    Returns
    -------
    d : datetime date
        date
    ampm : str, 'AM'|'PM'
        ampm indicator
    '''
    from datetime import date

    d = date(*t[:-1])     # current date, datetime date
    ampm = t[-1]          # current ampm, 'AM'|'PM'

    return d, ampm

def t2dt(t):
    ''' convert a t tuple to datetime
    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'

    Returns
    -------
    dt : datetime datetime
        datetime
    '''
    from datetime import datetime
    dlh = {'AM':0,'PM':12}

    dt = datetime(*t[:-1],dlh[t[-1]])

    return dt

def d2t(year,month,day,ampm):
    ''' convert year, month, day, ampm to a t tuple
    Parameters
    ----------
    year : int
        year
    month : int
        month
    day : int
        day
    ampm : str, 'AM'|'PM'
        ampm indicator

    Returns
    -------
    t : list, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    '''
    t = [year,month,day,ampm]
    return t

def dt2t(dt):
    ''' convert datetime to a t tuple
    Parameters
    ----------
    dt : datetime datetime
        datetime
    Returns
    -------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    '''
    dhl = {0:'AM',12:'PM'}

    t = [dt.year,dt.month,dt.day,dhl[dt.hour]]
    return t

def ftrange(firstday,lastday):
    ''' get datetime range for given first and last t tuples (both ends included)

    Parameters
    ----------
    firstday : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    lastday : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'

    Returns
    -------
    trange : pandas date range
        date range defined by firstday and lastday
    '''
    import pandas as pd

    trange = pd.date_range(t2dt(firstday),t2dt(lastday),freq='12h')
    return trange