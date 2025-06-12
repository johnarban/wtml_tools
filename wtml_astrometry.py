from astropy.io import fits
from io import BytesIO
from time import sleep
import requests
import json
from PIL import Image


def get_anet_session(api_key, anet=None, subid=None, jobid=None):
    new_anet = Astrometry(api_key)
    if anet is not None:
        new_anet.subid = anet.subid
        new_anet._jobid = anet._jobid
    else:
        if subid is not None:
            new_anet.subid = subid
        if jobid is not None:
            new_anet._jobid = jobid
    return new_anet


class Astrometry:
    api_key = None
    api_url = "http://nova.astrometry.net/api/"
    session = None
    subid = None
    _jobid = None
    status = None
    jobid_check_threshold = 10
    _hdu = None
    image_url = None
    _original_filename = None
    _image = None
    _image_id = None
    _calibration = None
    _status_url = None
    _results_url = None
    _user_id = None
    _post_data = None
    _post_response = None
    _user_image = None
    _user_image_url = None
    _header = None

    def __init__(self, api_key=None):
        if api_key is not None:
            self.api_key = api_key
        if self.api_key is None:
            raise ValueError("API key not set. Call with api_key=...")

        self.login()

    def login(self):
        self._post_data = {"request-json": json.dumps({"apikey": self.api_key})}
        R = requests.post(
            "http://nova.astrometry.net/api/login", data=self._post_data
        )
        self.session = R.json()["session"]
        self.login_session = R

    def submit(
        self,
        url,
        center_ra=None,
        center_dec=None,
        radius=None,
        scale_units="degwidth",
        scale_lower=0.1,
        scale_upper=20,
        block_till_done=False,
    ):
        upload_url = "http://nova.astrometry.net/api/url_upload"

        message = {
            "session": self.session,
            "allow_commercial_use": "n",
            "publicly_visible": "n",
            "parity": 2,
            "url": url,
            "scale_units": scale_units,
            "scale_lower": scale_lower,
            "scale_upper": scale_upper,
            "crpix_center": True,
            "tweak_order": 1,
        }

        self.image_url = url
        if center_ra is not None and center_dec is not None:
            radius = 10 or radius

        position = {
            "center_ra": center_ra,
            "center_dec": center_dec,
            "radius": radius,
        }
        if len([x for x in position.values() if x is not None]) == 3:
            message.update(position)

        R = requests.post(
            upload_url, data={"request-json": json.dumps(message)}
        )
        self._post_response = R
        if R.status_code != 200:
            raise ValueError(
                "Something went wrong with the upload. Status code: {} {}".format(
                    R.status_code, R.text
                )
            )

        self.subid = R.json()["subid"]
        print(self.status_url)

        if block_till_done:
            self.block_till_done()

    def block_till_done(self):
        print("submission id: {}".format(self.subid))
        # Print out the submission status
        print("Submission status:")
        for key, value in self.submission_status.items():
            print("\t{}: {}".format(key, value))
        if value is not None:
            # runs a loop until it has a job id
            if self.jobid is None:
                print(
                    "Waiting for jobid...(checks remaining: {} x 3sec)".format(
                        self.jobid_check_threshold
                    )
                )

            while (self.jobid is None) & (self.jobid_check_threshold >= 0):
                self.jobid_check_threshold -= 1
                sleep(3)

            if self.jobid_check_threshold <= 0:
                print("taking too long to get jobid")
                return

            while self.job_status is None:
                self.jobid
                sleep(1)

            while self.job_status == "processing":
                print("Job is processing")
                sleep(1)

            while self.job_status == "solving":
                print("Job is solving")
                sleep(1)

            print("Final Status: {}".format(self.job_status))

    @staticmethod
    def empty_list(l):
        return not bool([x for x in l if x is not None])

    def get_submission_status(self):
        submission_status_url = "http://nova.astrometry.net/api/submissions/{}"
        R = requests.get(submission_status_url.format(self.subid))

        return R

    @property
    def submission_status(self):
        R = self.get_submission_status()
        try:
            res = json.loads(R.text)
            self._user_id = R.json()["user"]
            return R.json()
        except:
            return None

    @property
    def status_url(self):
        if self.subid is not None:
            if self._status_url is None:
                self._status_url = (
                    "https://nova.astrometry.net/status/{}".format(self.subid)
                )
        return self._status_url

    @property
    def results_url(self):
        if self.subid is not None:
            if self._results_url is None:
                self._results_url = (
                    "https://nova.astrometry.net/user_images/{}".format(
                        self.user_image
                    )
                )
        return self._results_url

    @property
    def jobid(self):
        """caution this is recursive"""
        if self._jobid is None:
            sub = self.submission_status["jobs"]
            if self.empty_list(sub):
                return None
            else:
                self._jobid = sub[0]
                return self._jobid

        return self._jobid

    @jobid.setter
    def jobid(self, value):
        self._jobid = value

    @property
    def job_status(self):
        if self.jobid is not None:
            endoint = "http://nova.astrometry.net/api/jobs/{}"
            R = requests.get(endoint.format(self.jobid))
            return R.json()["status"]
        return None

    @property
    def image_id(self):
        if self._image_id is None:
            while self.empty_list(self.submission_status["images"]):
                sleep(1)
            self._image_id = self.submission_status["images"][0]
        return self._image_id

    @property
    def image(self):
        if self._image is None:
            api_image_url = "https://nova.astrometry.net/image/{}"
            R = requests.get(api_image_url.format(self.image_id))
            self._image = Image.open(BytesIO(R.content))
        return self._image

    def get_results(self):
        if self.jobid is not None:
            endoint = "http://nova.astrometry.net/api/jobs/{}/info"
            R = requests.get(endoint.format(self.jobid))
            return R.json()

    @property
    def calibration(self):
        if self.job_status is not None:
            endoint = "http://nova.astrometry.net/api/jobs/{}/calibration"
            R = requests.get(endoint.format(self.jobid))
            return R.json()

    def _wcs_file(self):
        if self.job_status is not None:
            endoint = "http://nova.astrometry.net/wcs_file/{}"
            R = requests.get(endoint.format(self.jobid))
            return R

    @property
    def header(self):
        if self.job_status is not None:
            if self._header is None:
                self._header = fits.Header.fromstring(self._wcs_file().text)
                return self._header
            else:
                return self._header

    def _new_fits_file(self):
        if self.job_status is not None:
            endoint = "http://nova.astrometry.net/new_fits_file/{}"
            R = requests.get(endoint.format(self.jobid))
            return R

    @property
    def new_fits(self):
        if self._hdu is None:
            if self.job_status is not None:
                bytesteam = BytesIO(self._new_fits_file().content)
                with fits.open(bytesteam) as hdul:
                    data = hdul[0].data
                    header = hdul[0].header
                    self._hdu = fits.PrimaryHDU(data=data, header=header)
        return self._hdu

    @property
    def hdu(self):
        if self._hdu is None:
            return self.new_fits
        else:
            return self._hdu

    @property
    def fits_filestream(self):
        return BytesIO(self.new_fits_file.content)

    def get_original_file_name(self):
        if self.image_url is not None:
            return self.image_url.split("/")[-1]
        else:
            if self.subid is not None:
                if self.jobid is None:
                    self.jobid
                result = self.get_results()
                if result is not None:
                    return result["original_filename"]

    @property
    def original_filename(self):
        if self._original_filename is None:
            self._original_filename = self.get_original_file_name()
        return self._original_filename

    @property
    def user_image(self):
        if self._user_image is None:
            user_image = self.submission_status["user_images"]
            if (user_image is not None) and (len(user_image) > 0):
                self._user_image = user_image[0]
                user_image_url = "https://nova.astrometry.net/user_images/{}"
                self._user_image_url = user_image_url.format(self._user_image)
        return self._user_image

    @property
    def user_image_url(self):
        if self._user_image_url is None:
            user_image = self.submission_status["user_images"]
            if (user_image is not None) and (len(user_image) > 0):
                self._user_image = user_image[0]
                user_image_url = "https://nova.astrometry.net/user_images/{}"
                self._user_image_url = user_image_url.format(self._user_image)
        return self._user_image_url
