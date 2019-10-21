|Build Status| |release| |Documentation Status| 
|HitCount| |MIT Licence| |Say Thanks!| |Donate|

Introduction
============

Nornir is a runtime support, that can be used to control performance and/or power consumption of parallel applications. The user can specify performance and power consumption requiremnets, and Nornir will enforce them by selecting the appropriate amount of resources to allocate to the application (e.g. number of cores, clock frequency, threads' mapping, etc...). The application will be monitored throughout its entire execution to provide such guarantees even in presence of workload fluctuations or phase changes. Nornir also allows the user to easily integrate custom decision policies inside the framework. For more information, please refer to the
`documentation <https://nornir.readthedocs.io/en/latest/index.html>`__.

Contributions
=============

Nornir has been mainly developed by Daniele De Sensi
(d.desensi.software@gmail.com).

The following people contributed to Nornir: 

- Daniele De Sensi (d.desensi.software@gmail.com): Main developer 
- Federico Umani: PredictorSMT predictor.

If you would like to contribute to Nornir development, for example by
adding new algorithms, please refer to the
`documentation <https://nornir.readthedocs.io/en/latest/contributing.html>`__.


.. |Build Status| image:: https://travis-ci.org/DanieleDeSensi/nornir.svg?branch=master
   :target: https://travis-ci.org/DanieleDeSensi/nornir
.. |release| image:: https://img.shields.io/github/release/danieledesensi/nornir.svg
   :target: https://github.com/danieledesensi/nornir/releases/latest
.. |Documentation Status| image:: https://readthedocs.org/projects/nornir/badge/?version=latest
   :target: https://nornir.readthedocs.io/en/latest/?badge=latest
.. |HitCount| image:: http://hits.dwyl.io/DanieleDeSensi/nornir.svg
   :target: http://hits.dwyl.io/DanieleDeSensi/nornir
.. |MIT Licence| image:: https://badges.frapsoft.com/os/mit/mit.svg?v=103
   :target: https://opensource.org/licenses/mit-license.php
.. |Say Thanks!| image:: https://img.shields.io/badge/Say%20Thanks-!-1EAEDB.svg
   :target: https://saythanks.io/to/DanieleDeSensi
.. |Donate| image:: https://img.shields.io/badge/Donate-PayPal-green.svg
   :target: http://paypal.me/DanieleDeSensi
.. _documentation: https://nornir.readthedocs.io/en/latest/