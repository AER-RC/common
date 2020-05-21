#!/usr/bin/env python

from __future__ import print_function

import os, sys
import requests

# part of the common_modules library
import utils

class apiZenodo:
  def __init__(self, inDict):
    """
    `zenodo_request.py -h`
    """

    self.url = inDict['url']
    self.tokenFile = inDict['access_token']
    utils.file_check(self.tokenFile)

    self.localDir = inDict['local_dir']
    self.title = inDict['dataset_title']
    self.sandbox = inDict['sandbox']

    # these two attributes are always the same for uploads and are
    # linked with the user account associated with the access token
    # https://developers.zenodo.org/?python#requests
    self.headers = {"Content-Type": "application/json"}
    self.depID = inDict['id_dep']

    # HTTP status codes for successful requests
    self.success = [200, 201, 202]

    # make sure we have data to upload
    utils.file_check(self.localDir)

    # upload URL will always be the same
    domain = 'sandbox.zenodo' if self.sandbox else 'zenodo'
    self.url = 'https://%s.org/api/deposit/depositions' % domain

    # create list of dictionaries for authors/institutions
    authors = inDict['creators']
    affils = inDict['affiliations']
    if len(authors) != len(affils):
      print('Number of creators and affiliations must be equal.')
      sys.exit(1)
    # endif authors

    # one dictionary per author
    self.creators = []
    for c, a in zip(authors, affils):
      self.creators.append({'name': c, 'affiliation': a})

    # all uploads will probably be datasets, but eventually we can
    # add more flexibility
    self.type = 'dataset'

    # communities can be changed in the web interface
    # for now, let's just throw it into the AER archive
    self.communities = [{'identifier': 'aer'}]
  # end constructor()

  def getKey(self):
    """
    Read in secret Zenodo API Access Token
    """

    import numpy as np
    self.token = str(np.loadtxt(self.tokenFile, dtype=str))
  # end getKey

  def zUpload(self):
    """
    Zenodo API upload request

    https://developers.zenodo.org/?python#quickstart-upload

    https://stackoverflow.com/a/22567429
    """

    import glob, json

    # list of files to upload
    upFiles = sorted(\
      glob.glob('%s/**/*.nc' % self.localDir, recursive=True))

    # create a new deposit -- for now, we'll do 1 per localDir
    # deposition starts off empty
    req = requests.post(self.url, headers=self.headers, \
      params={'access_token': self.token}, json={})
    self.depID = req.json()['id']

    for upFile in upFiles:
      base = os.path.basename(upFile)

      # name of file in Zenodo archive
      data = {'filename': base}

      # where Zenodo will get the data
      files = {'file': open(upFile, 'rb')}

      # upload file to dataset
      depURL = '%s/%s' % (self.url, self.depID)
      req = requests.post('%s/files' % depURL, \
        params={'access_token': self.token}, data=data, files=files)
      if req.status_code in self.success:
        print('Uploaded %s' % base)
      else:
        print('Status code for file upload: %d' % req.status_code)
        print(req.text)
      # end if success
    # end upload loop

    # attach associated metadata for entire dataset
    metadata = {'metadata': \
      {'title': self.title, 'upload_type': self.type, \
      'description': self.title, 'creators': self.creators, \
      'communities': self.communities}}
    req = requests.put(depURL, params={'access_token': self.token}, \
      data=json.dumps(metadata), headers=self.headers)
    if req.status_code in self.success:
      print('Metadata attached to deposit %d' % self.depID)
    else:
      print('Status code for metadata: %d' % req.status_code)
      print(req.text)
    # end if success

    # publish dataset
    req = requests.post('%s/actions/publish' % depURL, \
      params={'access_token': self.token})
    if req.status_code in self.success:
      print('Data published to deposit %d' % self.depID)
    else:
      print('Status code for publication: %d' % req.status_code)
      print(req.text)
    # end if success

  # end zUpload
# end apiZenodo

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(\
    description='Utilize the Zenodo API to upload ' + \
    'files from the service. Downloads are recommended to be ' + \
    'done with https://gitlab.com/dvolgyes/zenodo_get. Files are ' + \
    'are associated with the user as specified in access_token. ' + \
    'This is a reboot of my original script, which was designed ' + \
    'for NH3 retrievals.')
  parser.add_argument('access_token', type=str, \
    help='Text file with Zenodo Access Token in it.')
  parser.add_argument('--url', '-u', type=str, \
    default='https://zenodo.org/record/1979000/files/' + \
    'Austin_all_historical_data.csv?download=1', \
    help='URL of file to download')
  parser.add_argument('--upload', '-up', action='store_true', \
    help='Upload files in --local_dir to Zenodo.')
  parser.add_argument('--local_dir', '-d', type=str, default='2013', \
    help='Local directory of files to upload.')
  parser.add_argument('--dataset_title', '-t', type=str, \
    default='AER Line File Parameters', \
    help='Name for dataset. This will also be the dataset ' + \
    'description, so it should should be descriptive enough for ' + \
    'end user to understand what is is in it.')
  parser.add_argument('--creators', '-c', nargs='+', type=str, \
    default=['Karen Cady-Pereira', 'Matthew Alvardo', \
      'Mike Iacono', 'Eli Mlawer', 'Rick Pernak'], \
    help='List of each author responsible for generating dataset.')
  parser.add_argument('--affiliations', '-a', nargs='+', type=str, \
    default=['AER']*5, help='Institution for every creator.')
  parser.add_argument('--sandbox', '-s', action='store_true', \
    help='Upload to Zenodo Sandbox instead of publication space.')
  args = parser.parse_args()

  sys.exit('Work in progress!')

  zObj = apiZenodo(vars(args))
  zObj.getKey()

  zObj.zUpload()
# end main()
