import numpy as np
from datetime import datetime
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import CIRS, GCRS
from astropy.coordinates import SkyCoord, Distance
import p2api
from astroquery.simbad import Simbad
Simbad.add_votable_fields('pmra', 'pmdec', 'rv_value', 'sptype', 'otype')

now = datetime.now()
dt_string = now.strftime("%Y-%m-%dT%H:%M:%S")

all = ['check_item_exist', 'create_OB', 'ob_add_target', 'ob_add_description', 
       'p2api_SOFI', 'search_simbad', 'create_OB_telluric']


def check_item_exist(name, itemType, runContainerId, api):
    '''
    Check whether folder is in the container.
    
    Parameters
    ----------
    name : string
        The name of the item.
    itemType : string
        The type of the item, 'OB' or 'Folder'.
    runContainerId : int
        The ID of the container.
    api : p2api
    
    Returns
    -------
    it : dict 
        The information of the iterm if the item exists.  False, otherwise.
    '''
    items, _ = api.getItems(runContainerId)
    
    for it in items:
        if (it['itemType'] == itemType) & (it['name'] == name):
            return it
        
    return False


def create_OB(ob_name, folder_name, runContainerId, api):
    '''
    Create an OB or replace an OB in a folder.
    
    Parameters
    ----------
    ob_name : string
        The OB name.
    folder_name : string
        The folder name.
    runContainerId : int
        The ID of the container.
    api : p2api
    
    Returns
    -------
    ob : dict
        The OB information.
    '''
    folder_info = check_item_exist(folder_name, 'Folder', runContainerId, api)

    if folder_info:
        folderId = folder_info['containerId']
    else:
        folder, folderVersion = api.createFolder(runContainerId, folder_name)
        folderId = folder['containerId']
    
    ob_info = check_item_exist(ob_name, 'OB', folderId, api)
    if ob_info:
        obId = ob_info['obId']
        ob, obVersion = api.getOB(obId)
        api.deleteOB(obId, obVersion)
    ob, obVersion = api.createOB(folderId, ob_name)
    return ob, obVersion
    
    
def ob_add_target(ob, name, ra, dec, pma, pmd):
    '''
    Add the target information.
    
    Parameters
    ----------
    ob : dict
        The OB information.
    name : string
        The target name.
    ra : string
        The R.A. in hour angle, HH:MM:SS.SSS.
    dec : string
        The Decl. in degree, DD:MM:SS.SSS.
    pma : float
        The proper motion of R.A., in milliarcsec.
    pmd : float
        The proper motion of Decl., in milliarcsec.
    parallax : float
        The parallax, in milliarcsec.
    
    Returns
    -------
    ob : dict
        The OB information.
    '''
    ob['target']['name'] = name
    ob['target']['ra'] = ra
    ob['target']['dec'] = dec
    ob['target']['properMotionRa'] = np.round(pma * 1e-3, 6)
    ob['target']['properMotionDec'] = np.round(pmd * 1e-3, 6)
    return ob


def ob_add_description(ob, name, userComments):
    '''
    Parameters
    ----------
    ob : dict
        The OB information.
    name : string
        The observing description name.
    userComments : string
        The user comments.
    
    Returns
    -------
    ob : dict
        The OB information.
    '''
    ob['obsDescription']['name'] = name
    ob['obsDescription']['userComments'] = userComments
    return ob
    

class p2api_SOFI(object):
    '''
    Create OB on P2LS for SOFI.
    '''
    def __init__(self, prog_id, username, password, root_container_id=None, no_warning=False):
        '''
        Initiate the API.
        
        Parameters
        ----------
        prog_id : string 
            Program ID, e.g. "109.23CR.001".
        username : string
            User name.
        password : string
            Password.
        '''
        api = p2api.ApiConnection('production_lasilla', username, password)
        runList = api.getRuns()[0]
        pidList = [r['progId'] for r in runList]
        
        if prog_id in pidList:
            nidx = pidList.index(prog_id)
        else:
            raise ValueError('Cannot find {0} in {1}!'.format(prog_id, pidList))
        
        self.api = api
        self.prog_id = prog_id
        self._run = runList[nidx]
        self._runContainerId = runList[nidx]['containerId']
        self.rootDict = dict(Root=(runList[nidx], 0))
        if root_container_id is None:
            self._rootContainterId = self._runContainerId
        else:
            self._rootContainterId = _rootContainterId
        
        self.update_content()
        
        
    def add_SOFI_img_acq_MoveToPixel(self, name, folder_name=None, dit=1, ndit=1,
                                     filt1='Ks', filt2='open', imode='LARGE_FIELD_IMAGING', 
                                     offsetalpha=0, offsetdelta=0, addvelalpha=0, 
                                     addveldelta=0, offangle=0, offset=False, preset=True, 
                                     save=False):
        '''
        Add the SOFI_img_acq_MoveToPixel acquisition template.
        '''
        api = self.api
        folder, folderVersion = self.get_folder(folder_name)
        odict = folder.get('OBs', None)
        if odict is None:
            raise Exception('The folder ({}) does not contain any OB!'.format(folder_name))
        
        assert name in odict, 'Cannot find the OB in {}!'.format(folder_name)
        ob, obVersion = odict[name]
        
        obid = ob['obId']
        try:
            acqTpl, acqTplVersion = api.createTemplate(obid, 'SOFI_img_acq_MoveToPixel')
        except p2api.P2Error as e:
            raise p2api.P2Error(e)
        
        pars = acqTpl['parameters']
        pList = [dit, ndit, filt1, filt2, imode, offsetalpha, offsetdelta, 
                 addvelalpha, addveldelta, offangle, offset, preset, save]
        pdict = {}
        for loop, pv in enumerate(pList):
            pdict[pars[loop]['name']] = pv
        acqTpl, acqTplVersion  = api.setTemplateParams(obid, acqTpl, pdict, acqTplVersion)
        
        # update the OB version
        odict[name] = api.getOB(ob['obId'])
        
        
    def add_SOFI_img_acq_MoveToSlit(self, name, folder_name=None, dit=1, ndit=1, 
                                    filt1='Ks', filt2='open', whichslit='long_slit_1', 
                                    offsetalpha=0, offsetdelta=0, addvelalpha=0, 
                                    addveldelta=0, offangle=0, offset=False, preset=True, 
                                    save=False):
        '''
        Add the SOFI_img_acq_MoveToSlit acquisition template.
        '''
        api = self.api
        folder, folderVersion = self.get_folder(folder_name)
        odict = folder.get('OBs', None)
        if odict is None:
            raise Exception('The folder ({}) does not contain any OB!'.format(folder_name))
        
        assert name in odict, 'Cannot find the OB in {}!'.format(folder_name)
        ob, obVersion = odict[name]
        
        obid = ob['obId']
        try:
            acqTpl, acqTplVersion = api.createTemplate(obid, 'SOFI_img_acq_MoveToSlit')
        except p2api.P2Error as e:
            raise p2api.P2Error(e)
            
        pars = acqTpl['parameters']
        pList = [dit, ndit, filt1, filt2, whichslit, offsetalpha, offsetdelta, 
                 addvelalpha, addveldelta, offangle, offset, preset, save]
        pdict = {}
        for loop, pv in enumerate(pList):
            pdict[pars[loop]['name']] = pv
        acqTpl, acqTplVersion  = api.setTemplateParams(obid, acqTpl, pdict, acqTplVersion)
        
        # update the OB version
        odict[name] = api.getOB(ob['obId'])


    def add_SOFI_img_acq_Preset(self, name, folder_name=None, dit=1, ndit=1,
                                filt1='Ks', filt2='open', imode='LARGE_FIELD_IMAGING', 
                                addvelalpha=0, addveldelta=0, offangle=0):
        '''
        Add the SOFI_img_acq_Preset acquisition template.
        '''
        api = self.api
        folder, folderVersion = self.get_folder(folder_name)
        odict = folder.get('OBs', None)
        if odict is None:
            raise Exception('The folder ({}) does not contain any OB!'.format(folder_name))
        
        assert name in odict, 'Cannot find the OB in {}!'.format(folder_name)
        ob, obVersion = odict[name]
        
        obid = ob['obId']
        try:
            acqTpl, acqTplVersion = api.createTemplate(obid, 'SOFI_img_acq_Preset')
        except p2api.P2Error as e:
            raise p2api.P2Error(e)
            
        pars = acqTpl['parameters']
        pList = [dit, ndit, filt1, filt2, imode, addvelalpha, addveldelta, offangle]
        pdict = {}
        for loop, pv in enumerate(pList):
            pdict[pars[loop]['name']] = pv
        acqTpl, acqTplVersion  = api.setTemplateParams(obid, acqTpl, pdict, acqTplVersion)
        
        # update the OB version
        odict[name] = api.getOB(ob['obId'])
        
        
    def add_SOFI_img_obs_AutoJitter(self, name, folder_name=None, exp_name='SOFI', 
                                    dit=1, ndit=1, win_nx=1024, win_ny=1024, startx=1, 
                                    starty=1, nexpo=1, filt1='Ks', filt2='open', 
                                    imode='LARGE_FIELD_IMAGING', offset=False, 
                                    jitter_width=40, return_org=True):
        '''
        Add the SOFI_img_obs_AutoJitter exposure template.
        '''
        api = self.api
        folder, folderVersion = self.get_folder(folder_name)
        odict = folder.get('OBs', None)
        if odict is None:
            raise Exception('The folder ({}) does not contain any OB!'.format(folder_name))
        
        assert name in odict, 'Cannot find the OB in {}!'.format(folder_name)
        ob, obVersion = odict[name]
        
        obid = ob['obId']
        try:
            acqTpl, acqTplVersion = api.createTemplate(obid, 'SOFI_img_obs_AutoJitter')
        except p2api.P2Error as e:
            raise p2api.P2Error(e)
        
        pars = acqTpl['parameters']
        pList = [exp_name, dit, ndit, win_nx, win_ny, startx, starty, nexpo, 
                 filt1, filt2, imode, offset, jitter_width, return_org]
        pdict = {}
        for loop, pv in enumerate(pList):
            pdict[pars[loop]['name']] = pv
        acqTpl, acqTplVersion  = api.setTemplateParams(obid, acqTpl, pdict, acqTplVersion)
        
        # update the OB version
        odict[name] = api.getOB(ob['obId'])
        

    def add_SOFI_img_obs_AutoJitterRot(self, name, folder_name=None, exp_name='SOFI', 
                                       dit=1, ndit=1, win_nx=1024, win_ny=1024, startx=1, 
                                       starty=1, nexpo=1, filt1='Ks', filt2='open', 
                                       imode='LARGE_FIELD_IMAGING', offset=False, 
                                       jitter_width=40, return_org=True, offangle=[]):
        '''
        Add the SOFI_img_obs_AutoJitterRot exposure template.
        '''
        api = self.api
        folder, folderVersion = self.get_folder(folder_name)
        odict = folder.get('OBs', None)
        if odict is None:
            raise Exception('The folder ({}) does not contain any OB!'.format(folder_name))
        
        assert name in odict, 'Cannot find the OB in {}!'.format(folder_name)
        ob, obVersion = odict[name]
        
        obid = ob['obId']
        try:
            acqTpl, acqTplVersion = api.createTemplate(obid, 'SOFI_img_obs_AutoJitterRot')
        except p2api.P2Error as e:
            raise p2api.P2Error(e)
        
        pars = acqTpl['parameters']
        pList = [exp_name, dit, ndit, win_nx, win_ny, startx, starty, nexpo, 
                 filt1, filt2, imode, offset, jitter_width, return_org, offangle]
        pdict = {}
        for loop, pv in enumerate(pList):
            pdict[pars[loop]['name']] = pv
        acqTpl, acqTplVersion  = api.setTemplateParams(obid, acqTpl, pdict, acqTplVersion)
        
        # update the OB version
        odict[name] = api.getOB(ob['obId'])
        

    def add_SOFI_spec_obs_AutoNodNonDestr(self, name, folder_name=None, exp_name='SOFI', 
                                          dit=1, ndit=1, ndsamples=4, nsamppix=4, 
                                          win_nx=1024, win_ny=1024, startx=1, starty=1, 
                                          smode='LONG_SLIT_K', whichslit='long_slit_1', 
                                          offset=False, jitter_width=40, return_org=True, 
                                          nodthrow=60, nint=1, nabcycles=1):
        '''
        Add the SOFI_spec_obs_AutoNodNonDestr exposure template.
        '''
        api = self.api
        folder, folderVersion = self.get_folder(folder_name)
        odict = folder.get('OBs', None)
        if odict is None:
            raise Exception('The folder ({}) does not contain any OB!'.format(folder_name))
        
        assert name in odict, 'Cannot find the OB in {}!'.format(folder_name)
        ob, obVersion = odict[name]
        
        obid = ob['obId']
        try:
            acqTpl, acqTplVersion = api.createTemplate(obid, 'SOFI_spec_obs_AutoNodNonDestr')
        except p2api.P2Error as e:
            raise p2api.P2Error(e)
        
        pars = acqTpl['parameters']
        pList = [exp_name, dit, ndit, ndsamples, nsamppix, win_nx, win_ny, startx, 
                 starty, smode, whichslit, offset, jitter_width, return_org, nodthrow, 
                 nint, nabcycles]
        pdict = {}
        for loop, pv in enumerate(pList):
            pdict[pars[loop]['name']] = pv
        acqTpl, acqTplVersion  = api.setTemplateParams(obid, acqTpl, pdict, acqTplVersion)
        
        # update the OB version
        odict[name] = api.getOB(ob['obId'])
        
        
    def create_folder(self, name, container_id=None, overwrite=False):
        '''
        Create folder.
        
        Parameters
        ----------
        name : string
            Folder name.
        runContainerId (optional) : string
            The run container ID.  If not provided, create the folder directly under the run.
        overwrite : bool
            Overwrite the existing OB if True.
        '''
        api = self.api
        fdict = self.folderDict
        
        if container_id is None:
            container_id = self._rootContainterId
        
        if name in fdict:
            if overwrite:
                folder, folderVersion = fdict[name]
                self.delete_folder(folder['containerId'], force=True)
            else:
                raise Exception('The folder ({}) exists!'.format(name))
                
        fdict[name] = api.createFolder(container_id, name)
        return fdict[name]
    
        
    def create_OB(self, name, folder_name=None, overwrite=False):
        '''
        Create an OB.
        
        Parameters
        ----------
        name : string
            OB name.
        folder_name : string 
            Folder name.
        overwrite : bool
            Overwrite the existing OB if True.
        '''
        api = self.api
        try:
            folder, folderVersion = self.get_folder(folder_name)
        except AssertionError:
            folder, folderVersion = self.create_folder(folder_name)
        containerId = folder['containerId']
        
        odict = folder.get('OBs', None)
        if odict is None:
            folder['OBs'] = {}
            odict = folder['OBs']
            
        if name in odict:
            if overwrite:
                ob, obVersion = odict[name]
                api.deleteOB(ob['obId'], obVersion)
            else:
                raise Exception('The OB ({}) exists!'.format(name))
        
        odict[name] = api.createOB(containerId, name)
        return odict[name]
    
    
    def create_rootFolder(self, name, overwrite=False):
        '''
        Create a root folder and work on it.
        '''
        api = self.api
        containerId = self._runContainerId
        
        items = api.getItems(containerId)[0]
        for it in items:
            if (it['itemType'] == 'Folder') & (it['name'] == name):
                if overwrite:
                    self.delete_folder(it['containerId'], force=True)
                else:
                    raise Exception('The folder ({}) already exists!'.format(name))
        
        fd, fdV = api.createFolder(containerId, name)
        self.rootDict['Root'] = [fd, fdV]
        self.set_rootContainterId(fd['containerId'])
        return fd, fdV
    
    
    def delete_folder(self, folder_id, force=False):
        '''
        Delete the folder.
        
        Parameters 
        ----------
        folder_id : int 
            Folder ID.
        force : bool
            Empty the folder and delete it if True.
        '''
        api = self.api
        
        if force:
            # Delete all of the content first
            items = api.getItems(folder_id)[0]
            for it in items:
                if it['itemType'] == 'OB':
                    ob, obV = api.getOB(it['obId'])
                    api.deleteOB(ob['obId'], obV)
                elif it['itemType'] == 'Folder':
                    fd, fdV = api.getContainer(it['containerId'])
                    self.delete_folder(fd['containerId'], force=force)
                else:
                    raise KeyError('Cannot recognize this type ({})!'.format(it['itemType']))
            fd, fdV = api.getContainer(folder_id)
            api.deleteContainer(folder_id, fdV)
        else:
            fd, fdV = api.getContainer(folder_id)
            api.deleteContainer(folder_id, fdV)

    
    def get_folder(self, folder_name=None):
        '''
        Get the folder information.
        
        Parameters
        ----------
        folder_name (optional) : string 
            Folder name.
        '''
        if folder_name is None:
            folder_name = 'Root'
            folder, folderVersion = self.rootDict[folder_name]  # Root folder 
            #folderVersion = None
        else:
            assert folder_name in self.folderDict, 'Cannot find the foler ({})!'.format(folder_name)
            folder, folderVersion = self.folderDict[folder_name]
            
        return folder, folderVersion
        
        
    def save_OB(self, name, folder_name=None):
        '''
        Save an OB.
        
        Parameters
        ----------
        name : string
            OB name.
        folder_name : string 
            Folder name.
        '''
        api = self.api
        folder, folderVersion = self.get_folder(folder_name)
        odict = folder.get('OBs', None)
        if odict is None:
            raise Exception('The folder ({}) does not contain any OB!'.format(folder_name))
        
        assert name in odict, 'Cannot find the OB in {}!'.format(folder_name)
        ob, obVersion = odict[name]
            
        api.saveOB(ob, obVersion)
        odict[name] = api.getOB(ob['obId'])
        return odict[name]
    
    
    def set_rootFolder(self, name, no_warning=False):
        '''
        Set a root folder and work on it.
        '''
        api = self.api
        containerId = self._runContainerId
        
        items = api.getItems(containerId)[0]
        fd = None
        for it in items:
            if (it['itemType'] == 'Folder') & (it['name'] == name):
                fd, fdV = api.getContainer(it['containerId'])
                self.rootDict['Root'] = [fd, fdV]
        
        if fd is None:
            raise Exception('Cannot find a folder named {}!'.format(name))
        else:
            self.set_rootContainterId(fd['containerId'], no_warning)
        
    
    def set_rootContainterId(self, root_container_id, no_warning=False):
        '''
        Set the rootContainterId.
        
        Parameters
        ----------
        root_container_id : int 
            ID of the root container.
        '''
        self._rootContainterId = root_container_id
        self.update_content(no_warning)
    
    
    def update_content(self, no_warning=False):
        '''
        Update the content of this run.
        '''
        api = self.api
        
        # Put in the existing folders and OBs. I do not care the duplication at this point.
        self.folderDict = {}
        items = api.getItems(self._rootContainterId)[0]
        for it in items:
            name = it['name']
            if it['itemType'] == 'OB':
                if name in self.rootDict:
                    if not no_warning:
                        raise Warning('OB "{}" already exists in Root!'.format(name))
                self.rootDict[name] = api.getOB(it['obId'])
            elif it['itemType'] == 'Folder':
                if name in self.folderDict:
                    if not no_warning:
                        raise Warning('Folder "{}" already exists in Root!'.format(name))
                
                self.folderDict[name] = api.getContainer(it['containerId'])
                folder, folderVersion = self.folderDict[it['name']]
                folder['OBs'] = {}
                
                items_in_folder = api.getItems(folder['containerId'])[0]
                for it_in_folder in items_in_folder:
                    if it_in_folder['itemType'] == 'OB':
                        folder['OBs'][it_in_folder['name']] = api.getOB(it_in_folder['obId'])
            else:
                raise KeyError('Cannot recognize this type ({})!'.format(it['itemType']))


    def verifyOB(self, ob):
        '''
        Verify OB.
        '''
        self.api.verifyOB(ob['obId'], True)
        

def search_simbad(name):
    '''
    Find the target name from Simbad.
    
    Parameters
    ----------
    name : str
        Name of the target.
    
    Returns
    -------
    ra_hms, dec_hms  : str
        Coordinates, (HH:MM:SS, DD:MM:SS).
    pma, pmd : str
        Proper motion, units: arcsec / year
    '''
    result_table = Simbad.query_object(name)
    
    if result_table is None:
        raise Exception('The target ({}) is not found!'.format(name))
    
    ra = result_table['RA'][0]
    dec = result_table['DEC'][0]
    pma = result_table['PMRA'][0] * 1e-3
    pmd = result_table['PMDEC'][0] * 1e-3
    c = read_coordinate(ra, dec)
    ra_hms, dec_dms = get_coord_colon(c)
    return ra_hms, dec_dms, pma, pmd
    
    
def create_OB_telluric(name, sofi, folder_name='tmp', dit_acq=10, ndit_acq=1, dit=10, ndit=1, save=True):
    '''
    Create an OB for the telluric star.
    
    Parameters
    ----------
    name : str
        Name of the telluric star.
    sofi : p2api_SOFI
    folder_name : str (default: tmp)
        Name of the folder to save the OB.
    dit_acq : float (default: 10)
    ndit_acq : float (default: 1)
    dit : float (default: 10)
    ndit : float (default: 1)
    save : bool (default: True)
    '''
    ra, dec, pma, pmd = search_simbad(name)
    
    tell_name = name.replace(' ', '_')
    ob_name = 'telluric_{}'.format(tell_name)
    ob, obVersion = sofi.create_OB(ob_name, folder_name=folder_name, overwrite=True)
    target = ob['target']
    target['name'] = name
    target['ra'] = ra
    target['dec'] = dec
    target['properMotionRa'] = pma
    target['properMotionDec'] = pmd
    obsDes = ob['obsDescription']
    obsDes['name'] = name
    obsDes['userComments'] = 'Telluric star for {}'.format(folder_name)
    ob, obV = sofi.save_OB(ob_name, folder_name=folder_name)
    
    sofi.add_SOFI_img_acq_MoveToSlit(ob_name, folder_name=folder_name, dit=dit_acq, ndit=ndit_acq, save=save)
    sofi.add_SOFI_spec_obs_AutoNodNonDestr(ob_name, folder_name=folder_name, dit=dit, ndit=ndit, nint=1, nabcycles=1)
    sofi.verifyOB(ob)
    

def read_coordinate(ra, dec):
    '''
    Read in the coordinate, either in degree or hourangle. Only use ICRS frame.
    
    Parameters
    ----------
    ra : float or string
        The right ascension (degree or HH:MM:SS).
    dec : float or string
        The declination (degree or DD:MM:SS).
    
    Returns
    -------
    c : SkyCoord
        The coordinate object.
    '''
    if isinstance(ra, str):
        assert isinstance(dec, str)
        c = SkyCoord('{0} {1}'.format(ra, dec), frame='icrs', unit=(u.hourangle, u.deg))
    else:
        c = SkyCoord(ra, dec, frame='icrs', unit='deg')
    return c
    
    
def cal_coord_motion(c, pma=None, pmd=None, plx=None, radvel=None,
                     time_ref='J2000.0', time_cal='now', frame=GCRS):
    """
    Calculate the current coordinates considering the motion of the source.

    Parameters
    ----------
    c : Astropy SkyCoord
        The coordinate to calculate the motion.
    pma (optional) : Astropy Quantity (angular velocity)
        The proper motion of RA.
    pmd (optional) : Astropy Quantity (angular velocity)
        The proper motion of DEC.
    plx (optional) : Astropy Quantity (angle)
        The parallex.
    radvel (optional) : Astropy Quantity (velocity)
        The radial velocity.
    time_ref : string (default: 'J2000.0')
        The reference time of the input data.
    time_cal : string (default: 'now')
        The time to calculate the motion of the coordinate.
    frame : coordinate systems (default: GCRS)
        The coordinate system, which need to consider the orbit of the Earth
        around the Sun.

    Returns
    -------
    c_c : Astropy SkyCoord
        The calculated coordinate.
    """
    if pma is None:
        pma = 1e-16*u.mas/u.yr
    if pmd is None:
        pmd = 1e-16*u.mas/u.yr
    if plx is None:
        plx = 0 * u.mas
    if radvel is None:
        radvel = 0 * u.km/u.s
    c_m = SkyCoord(ra=c.ra, dec=c.dec, pm_ra_cosdec=pma, pm_dec=pmd,
                   distance=Distance(parallax=plx), radial_velocity=radvel,
                   frame=c.frame.name, obstime=time_ref)
    if time_cal == 'now':
        time_cal = Time.now()
    else:
        time_cal = Time(time_cal)
    c_c = c_m.apply_space_motion(time_cal).transform_to(GCRS(obstime=time_cal))
    return c_c


def cal_offset(c, c_ref):
    """
    Calculate the coordinate offset betweenn the two coordinates.

    Parameters
    ----------
    c : Astropy SkyCoord
        The coordinate to calculate the offset.
    c_ref : Astropy SkyCoord
        The reference coordinate.

    Returns
    -------
    (delta_ra, delta_dec) : (float, float)
        The offsets of right ascension and declination, units: arcsec.
    """
    sep = c_ref.separation(c)
    pa = c_ref.position_angle(c)
    delta_ra = (sep * np.sin(pa)).to('arcsec').value
    delta_dec = (sep * np.cos(pa)).to('arcsec').value
    return (delta_ra, delta_dec)
    
    
def get_coord_colon(c):
    '''
    Convert the coordinate to the (HH:MM:SS.SSS, DD:MM:SS.SSS) format.
    
    Parameters
    ----------
    c : SkyCoord
        The coordinate object.
    
    Returns
    -------
    ra_colon, dec_colon : string
        The coordinates in HH:MM:SS.SSS and DD:MM:SS.SSS
    '''
    try:
        ra, dec = c.to_string('hmsdms').split(' ')
    except ValueError:
        return '00:00:00.000', '00:00:00.000'
    ra_h = ra[:2]
    ra_m = ra[3:5]
    ra_s = float(ra[6:-1])
    dec_d, dec_tmp = dec.split('d')
    dec_m, dec_s = dec_tmp.split('m')
    dec_s = float(dec_s[:-1])
    ra_colon = '{0}:{1}:{2}'.format(ra_h, ra_m, '{0:.3f}'.format(ra_s).zfill(6))
    dec_colon = '{0}:{1}:{2}'.format(dec_d, dec_m, '{0:.3f}'.format(dec_s).zfill(6))
    return ra_colon, dec_colon
    
    
def spec_params(Kagn):
    if Kagn <= 12.5:  # 18 min
        dit = 30
        ndit = 3
        nint = 2
        nabc = 3
    elif (Kagn > 12.5) & (Kagn <= 13.5):  # 48 min
        dit = 30
        ndit = 3
        nint = 2
        nabc = 8
    elif (Kagn > 13.5) & (Kagn <= 16):  # 3.2 hour
        dit = 60
        ndit = 2
        nint = 2
        nabc = 15
    else:
        raise Exception('Too faint!!')
    return dit, ndit, nint, nabc

        
def imag_on_params(Kagn):
    '''
    Parameters for imaging OB for on-axis targets
    '''
    if (Kagn > 10) & (Kagn <= 11):
        dit = 2
        ndit = 2
        nexpo = 5
    elif (Kagn > 11) & (Kagn <= 12):
        dit = 4
        ndit = 2
        nexpo = 5
    elif (Kagn > 12) & (Kagn <= 13.5):
        dit = 8
        ndit = 2
        nexpo = 5
    else:
        raise Warning('Kagn ({}) Out of the useful range!'.format(Kagn))
    return dit, ndit, nexpo


def imag_off_params(Kagn):
    '''
    Parameters for imaging OB for off-axis targets
    '''
    if (Kagn <= 14):
        dit = 4
        ndit = 12
        nexpo = 5
    elif (Kagn > 14) & (Kagn <= 16):
        dit = 4
        ndit = 18
        nexpo = 5
    else:
        raise Warning('Kagn ({}) Out of the useful range!'.format(Kagn))
    return dit, ndit, nexpo