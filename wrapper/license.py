
import logging
import requests
from BeautifulSoup import BeautifulSoup as BS

logger = logging.getLogger(__name__)

class AvailableLicense(object):
    def __init__(self, url = 'https://caeliveusage.jpl.nasa.gov/live_usage_include.cfm'):
        """Get available license for a CAE software

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        # URL
        self.url = url

        # Logger
        self.logger = logging.getLogger(__name__)

        return

    def GetRequest(self):
        """ """
        # Read the get request
        self.logger.debug('Sending request to %s' %self.url)

        r = requests.get(self.url)

        try:
            soup = BS(r.content)
        except:
            raise ValueError('Unable to read the HTML. Stopping.')

        return soup

    def PostRequest(self, vendor, status, sort):
        """Get post request value """
        self.logger.debug('Posting request vendor : %s  status : %s  sort : %s'
                                                    %(vendor, status, sort))

        r = requests.post(self.url, data = {'vendor': vendor, 'status': status,
                                                            'sort': sort})
        try:
            soup = BS(r.content)
        except:
            raise ValueError('Unable to read the HTML. Stopping.')

        return soup

    def GetVendors(self):
        """Get list of CAE software vendors"""
        self.logger.debug('Get vendor list')

        soup = self.GetRequest()

        # Save vendor name and value in a dictionary
        vendors = {}
        ss = soup.find('select', attrs = {'name': 'vendor'})
        oo = ss.fetch('option')
        for item in oo[1:]:
            vendors[item.string.strip()] = item['value']

        return vendors

    def GetStatusValues(self):
        """Get status values """
        self.logger.debug('Get status values')

        soup = self.GetRequest()

        # Save status values and value in a dictionary
        statuses = {}
        ss = soup.find('select', attrs = {'name': 'status'})
        oo = ss.fetch('option')
        for item in oo[1:]:
            statuses[item.string] = item['value']

        return statuses

    def GetSortByValues(self):
        """Get sort by values """
        self.logger.debug('Get sort by values')

        soup = self.GetRequest()

        # Save status values and value in a dictionary
        sortbys = {}
        ss = soup.find('select', attrs = {'name': 'sort'})
        oo = ss.fetch('option')
        for item in oo[1:]:
            sortbys[item.string] = item['value']

        return sortbys

    def GetVendorSoftwares(self, vendor):
        """Get list of software for a particular vendor"""
        self.logger.debug('Get software for vendor : %s' %vendor)

        start_string = 'Users of'
        stop_string = ':'

        try:
            vendorval = self.GetVendors()[vendor]
        except KeyError:
            return

        self.logger.debug('Vendor value : %s' %vendorval)

        soup = self.PostRequest(vendorval, 'all', 'user')

        ss = soup.findAll('div', attrs = {'class': 'featurehead'})

        softwares = []
        for item in ss:
            softstr = item.string
            softwares.append(softstr[softstr.find(start_string) + 8:
                                     softstr.find(stop_string)])

        return softwares

    def NumLicenseAvail(self, vendor, software, status = 'all', sort = 'user'):
        """Determine number of available licenses for the software"""
        self.logger.debug('Get number of licenses available for %s' %software)

        start_string = 'Total of'
        stop_string = 'license'

        try:
            vendorval = self.GetVendors()[vendor]
        except KeyError:
            return

        self.logger.debug('Vendor value : %s' %vendorval)

        soup = self.PostRequest(vendorval, status, sort)

        ss = soup.findAll('div', attrs = {'class': 'featurehead'})

        (total_issued, total_inuse) = (None, None)
        for item in ss:
            if (software+':') in item.string:
                softstr = item.string
                total_issued = int(softstr[softstr.find(start_string) + 8:
                                                softstr.find(stop_string)])
                total_inuse = int(softstr[softstr.rfind(start_string) + 8:
                                                softstr.rfind(stop_string)])
                break

        return (total_issued, total_inuse, total_issued - total_inuse)

def CheckLicenseAvailbility(vendor, software, max_tries = 5):
    """Check software volume license and try max_tries before giving up.

    vendor : str
        Software vendor name

    software : str
        Software name

    max_tries : int
        Number of tries to acquire the license before giving up

    Returns
    -------
        : bool
        Whether license is acquired or not
    """
    ntries = 0

    al = AvailableLicense()
    while(ntries < max_tries):
        (issued, in_use, available) = al.NumLicenseAvail(vendor, software)

        if available == 0:
            logger.info('No license available for %s by %s. Trying again in 1 minute. Try #%d' %(software, vendor, ntries+1))
            sleep(60)
            ntries += 1
        else:
            logger.info('%d available licenses out of %d' %(available, issued))
            return True

    if ntries == max_tries:
        return False


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG, format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)

    al = AvailableLicense()
    #print(al.GetStatusValues())
    #print(al.GetSortByValues())
    #print(al.GetVendorSoftwares('Mathworks'))
    # print(al.GetVendors())
    # print(al.GetStatusValues())
    # print(al.GetSortByValues())
    #print(al.NumLicenseAvail('Sigmadyne', 'SIGFIT_STANDARD'))
    #print(al.NumLicenseAvail('Cullimore', 'ThermalDesktop'))
    print(al.NumLicenseAvail('MSC', 'CAMPUS'))
    print(al.NumLicenseAvail('Mathworks', ' MATLAB'))
    # print(al.GetVendorSoftwares('Sigmadyne'))
    # data = al.NumLicenseAvail('Cullimore', 'ThermalDesktop')

